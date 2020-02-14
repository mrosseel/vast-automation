from typing import Tuple

import toml

import main_vast
import do_compstars
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
import tqdm
import numpy as np
import do_aavso_report
import do_calibration
import reading
import utils
import argparse
import logging
import gc
import matplotlib as mplotlib
from collections import namedtuple
from multiprocessing import cpu_count, Manager
from comparison_stars import ComparisonStars
from functools import partial
from gatspy.periodic import LombScargleFast
from gatspy.periodic import TrendedLombScargle
from reading import trash_and_recreate_dir
from timeit import default_timer as timer
from star_description import StarDescription
from pathlib import Path
import pandas as pd
from pandas import DataFrame, Series

from star_metadata import SiteData, CatalogData, CompStarData

gc.enable()
mplotlib.use('Agg')  # needs no X server
Period = namedtuple('period', 'period origin')
TITLE_PAD = 40


def set_seaborn_style():
    sns.set_context("notebook", font_scale=4)
    sns.set_style("ticks")


def plot_lightcurve_raw(star: StarDescription, curve: DataFrame, chartsdir):
    logging.debug(f"Plotting raw lightcurve for {star.local_id}")
    return _plot_lightcurve(star, curve, chartsdir)


def plot_lightcurve_pa(star: StarDescription, curve: DataFrame, chartsdir, period: Period):
    logging.debug(f"Plotting phase adjusted lightcurve for {star.local_id}")
    convert_func = partial(phase_lock_lightcurve, period=period)
    return _plot_lightcurve(star, curve, chartsdir, f"_lightpa", convert_func, xlabel='phase-adjusted JD')


def _plot_lightcurve(star: StarDescription, curve: DataFrame, chartsdir, suffix=f"_light",
                     jd_adjusting_func=None, xlabel='JD'):
    try:
        star_id = star.local_id
        logging.debug(f"Plotting lightcurve for {star_id}")
        vsx_name, separation, extradata, filename_no_ext = utils.get_star_or_catalog_name(star, suffix=suffix)
        var_type = f"Type: {extradata['Type']}" if extradata is not None and 'Type' in extradata else ""
        save_location = Path(chartsdir, filename_no_ext + '.png')
        start = timer()
        upsilon_match = star.get_metadata('UPSILON') if star.has_metadata('UPSILON') else None
        upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
        end = timer()
        logging.debug(f"timing upsilon stuff {end - start}")
        coord = star.coords
        names = utils.get_star_names(star)
        catalog_title = f"{names[0]}" if names is not None and names is not star.local_id else ''
        plot_title = f"{catalog_title}\nStar {star.local_id}"

        if curve is None:
            logging.info(f"Curve is None for star {star_id}")
            return

        if jd_adjusting_func is None:
            curve['realJD'] = curve['floatJD']
        else:
            curve['realJD'] = jd_adjusting_func(curve['floatJD'])
        fig = plt.figure(figsize=(20, 12), dpi=150, facecolor='w', edgecolor='k')
        plt.errorbar(curve['realJD'], curve['realV'], yerr=curve['realErr'], linestyle='none', marker='o',
                     ecolor='gray',
                     elinewidth=1)

        plt.xlabel(xlabel, labelpad=TITLE_PAD)
        curve_max = curve['realV'].max()
        curve_min = curve['realV'].min()
        plot_max = curve_max + 0.1
        plot_min = curve_min - 0.1
        plt.ylim(plot_min, plot_max)
        nop = lambda *a, **k: None
        print("limits are here:", star_id, "min:", curve['realJD'].min(), "max", curve['realJD'].max(),
              "first 10", curve['realJD'][:10], "first 10 orig:", curve['floatJD'][:10], "len", len(curve['realJD']),
              "describe", curve['realJD'].describe()) \
            if curve['realJD'].max() == 0 else nop()

        plt.xlim(curve['realJD'].min(), curve['realJD'].max())
        plt.gca().invert_yaxis()
        # plt.ticklabel_format(style='plain', axis='x')
        plt.title(plot_title, pad=TITLE_PAD)
        plt.tight_layout()
        start = timer()
        fig.savefig(save_location)
        end = timer()
        logging.debug(f"timing saving fig {end - start}")
        plt.close(fig)
        plt.clf()
        return save_location
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        import traceback
        print(traceback.print_exc())
        logging.error(message)
        logging.error(f"Error during plot lightcurve: {star.local_id}")


def phase_lock_lightcurve(series: Series, period: Period):
    jds = np.sort(series.to_numpy())
    jds_norm_orig = np.subtract(jds, jds.min())
    jds_norm = np.diff(jds_norm_orig)
    jds_norm = np.mod(jds_norm, period.period)
    jds_norm = np.concatenate(([jds_norm_orig[0]], jds_norm))
    jds_norm = np.cumsum(jds_norm)
    return jds_norm


def plot_phase_diagram(star: StarDescription, curve: DataFrame, fullphasedir, suffix='', period: Period = None,
                       epoch: str = None, write_plot=True, filter_func=None):
    assert period is not None
    try:
        logging.debug(f"Starting plot phase diagram with {star} and {fullphasedir}")
        vsx_name, separation, extradata, filename_no_ext = utils.get_star_or_catalog_name(star,
                                                                                          suffix=f"_phase{suffix}")
        names = utils.get_star_names(star)
        catalog_title = f"{names[0]}" if names is not None and names is not star.local_id else ''

        save_location = Path(fullphasedir, filename_no_ext + '.png')
        upsilon_match = star.get_metadata('UPSILON')
        upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
        # print("Calculating phase diagram for", star)
        if curve is None:
            logging.info("Curve of star {} is None".format(star.local_id))
            return
        t_np = curve['floatJD']
        y_np = curve['realV'].to_numpy()
        dy_np = curve['realErr'].to_numpy()
        t_np_zeroed = shift_to_epoch(epoch, t_np)
        phased_t = np.mod(t_np_zeroed / period.period, 1)
        phased_lc = y_np[:]

        if filter_func is not None:
            phased_t, phased_lc = filter_func(phased_t, phased_lc)
        phased_t_final = np.append(phased_t.subtract(1), phased_t)
        phased_lc_final = np.append(phased_lc, phased_lc)
        phased_err = np.clip(np.append(dy_np, dy_np), -0.5, 0.5)  # error values are clipped to +0.5 and -0.5
        plt_result = _plot_phase_diagram(phased_t_final, phased_lc_final, phased_err, write_plot, save_location, star,
                                         catalog_title, period, upsilon_text)
        if not write_plot:
            return plt_result, t_np, y_np
        return save_location
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        import traceback
        print(traceback.print_exc())
        logging.error(message)
        logging.error(f"Error during plot phase: {star.local_id}")


def shift_to_epoch(epoch: float, t_np):
    """ shift the """
    if not epoch:
        return t_np
    assert isinstance(epoch, float)
    t_epoch_location = (np.abs(t_np - epoch)).idxmin()
    t_np_zeroed = t_np - t_np[t_epoch_location]
    return t_np_zeroed


def calculate_min_max_epochs(t_np, y_np):
    ymin_arg, ymax_arg = np.argmin(np.array(y_np)), np.argmax(np.array(y_np))
    epoch_min, epoch_max = t_np.iloc[ymin_arg], t_np.iloc[ymax_arg]
    ymin, ymax = y_np.iloc[ymin_arg], y_np.iloc[ymax_arg]
    return ymin, ymax, epoch_min, epoch_max


def write_toml(filename_no_ext, fullphasedir, period, star, points_removed, ymin, ymax, epoch_min, epoch_max):
    tomldict = {}
    tomldict['period'] = float(period.period)
    tomldict['period_origin'] = period.origin
    tomldict['range'] = f'{ymin:.1f}-{ymax:.1f} (LS)'
    tomldict['coords'] = [star.coords.ra.deg, star.coords.dec.deg]
    tomldict['points_removed'] = points_removed
    tomldict['our_name'] = f"{star.local_id}"
    if star.has_metadata('SITE'):
        sitedata: SiteData = star.get_metadata('SITE')
        assert sitedata is not None
        tomldict['var_type'] = sitedata.var_type
        tomldict['our_name'] = utils.get_star_names(star)
        tomldict['minmax'] = sitedata.minmax
        tomldict['min'] = sitedata.var_min
        tomldict['max'] = sitedata.var_max
        if sitedata.vsx_var_flag is not None:
            tomldict['vsx_var_flag'] = int(sitedata.vsx_var_flag)
        if sitedata.separation is not None:
            tomldict['separation'] = float(sitedata.separation)
        if sitedata.var_min and sitedata.var_max:
            tomldict['range'] = f'{sitedata.var_min:.2f}-{sitedata.var_max:.2f} ({sitedata.source})'
        if sitedata.period_err is not None:
            tomldict['period_err'] = sitedata.period_err
        if sitedata.epoch:
            tomldict['epoch'] = sitedata.epoch

    outputfile = f"{fullphasedir}/txt/{filename_no_ext}.txt"
    logging.debug(f"Writing toml to {outputfile}")
    toml.dump(tomldict, open(outputfile, "w"))


# plotting of 'double' phase diagram from -1 to 1
def _plot_phase_diagram(phased_t_final, phased_lc_final, phased_err, write_plot, save_location, star, catalog_title,
                        period: Period, upsilon_text):
    fig = plt.figure(figsize=(18, 16), dpi=80, facecolor='w', edgecolor='k')
    plt.xlabel("Phase", labelpad=TITLE_PAD)
    plt.ylabel("Magnitude", labelpad=TITLE_PAD)
    plt.title(
        f"{catalog_title}\nStar {star.local_id}, p: {period.period:.5f} d ({period.origin}) "
        # f"{upsilon_text}\n{utils.get_hms_dms_matplotlib(star.coords)}", pad=TITLE_PAD)
        f"{upsilon_text}", pad=TITLE_PAD)
    plt.ticklabel_format(style='plain', axis='x')
    plt.tight_layout()

    plt.gca().invert_yaxis()
    plt.errorbar(phased_t_final, phased_lc_final, yerr=phased_err, linestyle='none', marker='o', ecolor='gray',
                 elinewidth=1)
    if write_plot:
        logging.debug(f"Saving phase plot to {save_location}")
        fig.savefig(save_location, format='png')
        plt.close(fig)
        plt.clf()
    else:
        return plt


def write_compstars(star, filename_no_ext, fullphasedir, compstars, check_star):
    if compstars:
        outputfile = f"{fullphasedir}/txt/{filename_no_ext}_comps.txt"
        extra_id_ucac4 = utils.get_ucac4_of_sd(check_star.star_descriptions[0])
        with open(outputfile, "wt") as fp:
            fp.write(f"# K star: {extra_id_ucac4},{check_star.comp_catalogmags[0]:.3f},"
                     f"{check_star.comp_catalogerr[0]:.5f}\n")
            for sd, mag, err in zip(compstars.star_descriptions, compstars.comp_catalogmags, compstars.comp_catalogerr):
                ucac4_id = utils.get_ucac4_of_sd(sd)
                fp.write(f"{ucac4_id},{mag:.3f},{err:.5f}\n")


def read_vast_lightcurves(star: StarDescription, compstarproxy, star_result_dict, do_light, do_light_raw, do_phase,
                          do_aavso, aavso_limit, basedir: str, chartsdir: Path, phasedir: Path, aavsodir: Path):
    start = timer()
    if star.local_id not in star_result_dict:
        star_result_dict[star.local_id] = {}
    temp_dict = star_result_dict[star.local_id]
    if star.path is '':
        logging.debug(f"Path for {star.local_id} is empty")
        return
    if not do_light and not do_phase:
        logging.debug("Nothing to do, no charts or phase needed")

    logging.debug(
        f"Reading lightcurves for star {star.local_id} at path {star.path} for {star}...")
    # comp_mags = [x.vmag for x in comparison_stars]

    try:
        df = reading.read_lightcurve_vast(star.path)
        if df is None or len(df) == 0:
            logging.info(f"No lightcurve found for {star.path}")
            return
        comp_stars = compstarproxy.value
        filtered_compstars, check_star = do_compstars.filter_comparison_stars(star, comp_stars)
        comp_stars = None
        df['realV'], df['realErr'] = do_compstars.calculate_ensemble_photometry(
            df, filtered_compstars, do_compstars.weighted_value_ensemble_method)
        df['floatJD'] = df['JD'].astype(np.float)
        _, _, _, filename_no_ext = utils.get_star_or_catalog_name(star, suffix="")
        period, epoch = determine_period_and_epoch(df, star)
        df, points_removed = phase_dependent_outlier_removal(df, period)
        write_compstars(star, filename_no_ext, phasedir, filtered_compstars, check_star)
        write_toml(filename_no_ext, phasedir, period, star, points_removed,
                   *calculate_min_max_epochs(df['floatJD'], df['realV']))

        if do_phase and 'phase' not in star.result:
            temp_dict['phase'] = plot_phase_diagram(star, df.copy(), phasedir, period=period, epoch=epoch, suffix="")
        if do_light and 'lightpa' not in star.result:
            temp_dict['lightpa'] = plot_lightcurve_pa(star, df.copy(), chartsdir, period)
        if do_light_raw and 'light' not in star.result:
            temp_dict['light'] = plot_lightcurve_raw(star, df.copy(), chartsdir)
        if do_aavso and 'aavso' not in star.result:
            settings = toml.load('settings.txt')
            temp_dict['aavso'] = do_aavso_report.report(star, df.copy(), filtered_compstars, check_star, target_dir=aavsodir,
                        sitelat=settings['sitelat'], sitelong=settings['sitelong'], sitealt=settings['sitealt'],
                        camera_filter='V', observer=settings['observer'], chunk_size=aavso_limit)
        filtered_compstars = None
        star_result_dict[star.local_id] = temp_dict
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        logging.error(message)
        import traceback
        print(traceback.print_exc())
        logging.error(f"Exception during read_lightcurve for {star.path}, size JD: {len(df['floatJD'])},"
                      f"size V:  {len(df['realV'])}")

    end = timer()
    logging.debug(f"Full lightcurve/phase: {end - start}")


def phase_dependent_outlier_removal(df: DataFrame, period: Period) -> Tuple[DataFrame, int]:
    phased_t = np.fmod(df['floatJD'] / period.period, 1)
    # array of times rounded to 1 decimal, results in 11 buckens which cover the phase diagram from 0.0 to 1.0
    grouper = np.round(phased_t, 1)
    df_v_grouped = df.groupby(grouper)
    maskresult = pd.DataFrame()
    for bucket_label, bucket_all in df_v_grouped:
        bucket = bucket_all['realV']
        bucketmean = bucket.median()
        bucketstd = bucket.std()
        mask = abs(bucket - bucketmean) < 5 * bucketstd
        maskresult = maskresult.append(bucket_all[mask])
    return maskresult, len(df) - len(maskresult)


def determine_period_and_epoch(df: DataFrame, star: StarDescription) -> Tuple[Period, str]:
    if not star.has_metadata("SITE") or star.get_metadata("SITE").period is None:
        period: Period = calculate_ls_period_from_df(df.copy())
        logging.debug(f"Using LS period for star {star.local_id}: {period}")
        epoch = None
    else:
        sitedata = star.get_metadata("SITE")
        source = sitedata.source
        period: Period = Period(sitedata.period, source)
        epoch = sitedata.epoch
        logging.debug(f"Using {source} period for star {star.local_id}: {period}, epoch: {epoch}")
    return period, epoch


def calculate_ls_period_from_df(df: DataFrame) -> Period:
    return calculate_ls_period(df['floatJD'], df['realV'].to_numpy(), df['realErr'].to_numpy())


def calculate_ls_period(t_np, y_np, dy_np) -> Period:
    period_max = np.max(t_np) - np.min(t_np)
    if period_max <= 0.01:
        return
    ls = LombScargleFast(optimizer_kwds={'quiet': True, 'period_range': (0.01, period_max)},
                         silence_warnings=True, fit_period=True).fit(t_np, y_np, dy_np)
    period = ls.best_period
    return Period(period, "LS")
    # TODO test this trended lombscargle !
    # tmodel = TrendedLombScargle(optimizer_kwds={'quiet': True, 'period_range': (0.01, period_max)},
    #                             silence_warnings=True).fit(t_np, y_np)
    # period = tmodel.best_period


# reads lightcurves and passes them to lightcurve plot or phase plot
def run(star_descriptions, comp_stars: ComparisonStars, basedir: str, resultdir: str, phasepart: str, chartspart: str,
        aavso_part: str, do_light=False, do_light_raw=False, do_phase=True, do_aavso=False, aavsolimit=None,
        nr_threads=cpu_count(), desc="Writing light curve charts/phase diagrams"):
    chunk: int = 1  # max(1, len(star_descriptions) // nr_threads*10)
    set_seaborn_style()
    pool = mp.Pool(nr_threads, maxtasksperchild=10)
    phasedir = Path(resultdir, phasepart)
    chartsdir = Path(resultdir, chartspart)
    aavsodir = Path(resultdir, aavso_part)
    logging.debug(f"Using {nr_threads} threads for lightcurve, phase plotting and aavso reporting.")

    if do_phase:
        trash_and_recreate_dir(phasedir)
        trash_and_recreate_dir(Path(phasedir, Path('txt')))
    if do_light or do_light_raw:
        trash_and_recreate_dir(chartsdir)
    if do_aavso:
        trash_and_recreate_dir(aavsodir)
    with Manager() as manager:
        comp_stars_proxy = manager.Value('comp_stars', comp_stars)
        star_result_dict = manager.dict({})

        func = partial(read_vast_lightcurves, basedir=basedir, compstarproxy=comp_stars_proxy,
                       star_result_dict=star_result_dict, do_light=do_light, do_light_raw=do_light_raw,
                       do_phase=do_phase, do_aavso=do_aavso, aavso_limit=aavsolimit, phasedir=phasedir,
                       chartsdir=chartsdir, aavsodir=aavsodir)
        with tqdm.tqdm(total=len(star_descriptions), desc=desc, unit='stars') as pbar:
            for _ in pool.imap_unordered(func, star_descriptions, chunksize=chunk):
                pbar.update(1)
                pass
        pool.close()
        pool.join()
        stardict = main_vast.get_star_description_cache(star_descriptions)
        for key, value in star_result_dict.items():
            star_result = stardict[key].result
            for key2, value2 in value.items():
                if key2 not in star_result:
                    star_result[key2] = value2



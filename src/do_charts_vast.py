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
from pathlib import PurePath
from pandas import DataFrame, Series

from star_metadata import StarFileData

gc.enable()
mplotlib.use('Agg')  # needs no X server
Period = namedtuple('Period', 'period origin')
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
        vsx_title = '' if vsx_name is None else f"{vsx_name} Type: {extradata['Type']}"
        vsx_dist = '' if separation is None else f"(+/- {separation:.3f} deg)"
        save_location = PurePath(chartsdir, filename_no_ext + '.png')
        start = timer()
        upsilon_match = star.get_metadata('UPSILON') if star.has_metadata('UPSILON') else None
        upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
        end = timer()
        logging.debug(f"timing upsilon stuff {end - start}")
        coord = star.coords
        plot_title = f"{vsx_title}\nStar {star.local_id}  {vsx_dist}\n" \
                     f"{utils.get_hms_dms_matplotlib(coord)}{upsilon_text}"

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
        # if curve['realJD'].isna().any():
        #     print("we have a nan", jd_adjusting_func)
        #     print("star: ", star)
        #     if star.has_metadata("VSX"):
        #         print("extradata", star.get_metadata("VSX").extradata)
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
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        import traceback
        print(traceback.print_exc())
        logging.error(message)
        logging.error(f"Error during plot lightcurve: {star.local_id}")


def phase_lock_lightcurve(series: Series, period: Period):
    jds = series.to_numpy()
    jds_norm_orig = np.subtract(jds, jds.min())
    jds_norm = np.diff(jds_norm_orig)
    jds_norm = np.mod(jds_norm, period.period)
    jds_norm = np.concatenate(([jds_norm_orig[0]], jds_norm))
    jds_norm = np.cumsum(jds_norm)
    return jds_norm


def plot_phase_diagram(star: StarDescription, curve: DataFrame, fullphasedir, suffix='', period: Period = None,
                       write_plot=True, filter_func=None):
    assert period is not None
    try:
        logging.debug(f"Starting plot phase diagram with {star} and {fullphasedir}")
        vsx_name, _, extradata, filename_no_ext = utils.get_star_or_catalog_name(star, suffix=f"_phase{suffix}")
        catalog_title = f"{vsx_name}" if vsx_name is not None else ''

        save_location = PurePath(fullphasedir, filename_no_ext + '.png')
        upsilon_match = star.get_metadata('UPSILON')
        upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
        # print("Calculating phase diagram for", star)
        if curve is None:
            logging.info("Curve of star {} is None".format(star.local_id))
            return
        t_np = curve['floatJD']
        y_np = curve['realV'].to_numpy()
        dy_np = curve['realErr'].to_numpy()
        phased_t = np.fmod(t_np / period.period, 1)
        phased_lc = y_np[:]
        if filter_func is not None:
            phased_t, phased_lc = filter_func(phased_t, phased_lc)
        phased_t_final = np.append(phased_t.subtract(1), phased_t)
        phased_lc_final = np.append(phased_lc, phased_lc)
        phased_err = np.clip(np.append(dy_np, dy_np), -0.5, 0.5)  # error values are clipped to +0.5 and -0.5
        plt = _plot_phase_diagram(phased_t_final, phased_lc_final, phased_err, write_plot, save_location, star,
                                  catalog_title, period, upsilon_text)
        if not write_plot:
            return plt, t_np, y_np
        # write TXT file
        tomldict = {}
        ymin_arg, ymax_arg = np.argmin(y_np), np.argmax(y_np)
        epoch_min, epoch_max = t_np.iloc[ymin_arg], t_np.iloc[ymax_arg]
        ymin, ymax = y_np[ymin_arg], y_np[ymax_arg]
        var_range = f'{ymin:.1f}-{ymax:.1f}'
        epoch = None
        minmax = None
        if star.has_metadata('STARFILE'):
            starfiledata: StarFileData = star.get_metadata('STARFILE')
            assert starfiledata is not None
            if starfiledata.var_type is not None \
                and ("RR" in starfiledata.var_type or "W Uma" in starfiledata.var_type):
                if "RR" in starfiledata.var_type:
                    epoch = epoch_max
                    minmax = f"max: {ymax:.1f}"
                elif "W Uma" in starfiledata.var_type:
                    epoch = epoch_min
                    minmax = f"min: {ymin:.1f}"
                starfiledata.epoch = epoch
                starfiledata.minmax = minmax
            starfiledata.var_min = ymin
            starfiledata.var_max = ymax
            if starfiledata.period_err is not None:
                tomldict['period_err'] = float(starfiledata.period_err)
        tomldict['period'] = float(period.period)
        tomldict['period_origin'] = period.origin
        tomldict['min'] = float(ymin)
        tomldict['max'] = float(ymax)
        tomldict['range'] = var_range
        tomldict['coords'] = [star.coords.ra.deg, star.coords.dec.deg]
        if minmax:
            tomldict['minmax'] = minmax
        if epoch:
            tomldict['epoch'] = float(epoch)
        outputfile = f"{fullphasedir}/txt/{filename_no_ext}.txt"
        logging.debug(f"Writing toml to {outputfile}")
        toml.dump(tomldict, open(outputfile, "w"))
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        import traceback
        print(traceback.print_exc())
        logging.error(message)
        logging.error(f"Error during plot phase: {star.local_id}")


def _plot_phase_diagram(phased_t_final, phased_lc_final, phased_err, write_plot, save_location, star, vsx_title,
                        period: Period, upsilon_text):
    fig = plt.figure(figsize=(18, 16), dpi=80, facecolor='w', edgecolor='k')
    plt.xlabel("Phase", labelpad=TITLE_PAD)
    plt.ylabel("Magnitude", labelpad=TITLE_PAD)
    plt.title(
        f"{vsx_title}\nStar {star.local_id}, p: {period.period:.5f} d ({period.origin}) "
        f"{upsilon_text}\n{utils.get_hms_dms_matplotlib(star.coords)}", pad=TITLE_PAD)
    plt.ticklabel_format(style='plain', axis='x')
    # plt.title(f"Star {star} - {period.period}", pad=TITLE_PAD)
    plt.tight_layout()
    # plotting + calculation of 'double' phase diagram from -1 to 1
    plt.gca().invert_yaxis()
    plt.errorbar(phased_t_final, phased_lc_final, yerr=phased_err, linestyle='none', marker='o', ecolor='gray',
                 elinewidth=1)
    logging.debug(f"Saving phase plot to {save_location}")
    if write_plot:
        fig.savefig(save_location, format='png')
        plt.close(fig)
        plt.clf()
    else:
        return plt


def format_date(x, pos=None):
    thisind = np.clip(int(x + 0.5), 0, N - 1)
    return r.date[thisind].strftime('%Y-%m-%d')


def read_vast_lightcurves(star: StarDescription, compstarproxy, do_light, do_light_raw, do_phase, do_aavso,
                          aavso_limit, basedir: str, chartsdir: PurePath, phasedir: PurePath, aavsodir: PurePath):
    start = timer()
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
        filtered_compstars = do_compstars.get_star_compstars_from_catalog(star, comp_stars)
        comp_stars = None
        df['realV'], df['realErr'] = do_compstars.calculate_ensemble_photometry(
            df, filtered_compstars, do_compstars.weighted_value_ensemble_method)
        df['floatJD'] = df['JD'].astype(np.float)

        # remove errors
        # df = utils.reject_outliers_iqr(df, 'realV', 20)
        # logging.info(f"Rejected {old_size-len(df)} observations with iqr.")

        # calculate period for phase or light
        do_period = do_phase or do_light or do_light_raw
        if do_period:
            period = determine_period(df, star)
        if do_phase:
            plot_phase_diagram(star, df.copy(), phasedir, period=period, suffix="")
        if do_light:
            plot_lightcurve_pa(star, df.copy(), chartsdir, period)
        if do_light_raw:
            plot_lightcurve_raw(star, df.copy(), chartsdir)
        if do_aavso:
            settings = toml.load('settings.txt')
            do_aavso_report.report(star, df.copy(), filtered_compstars, target_dir=aavsodir,
                                   sitelat=settings['sitelat'], sitelong=settings['sitelong'],
                                   sitealt=settings['sitealt'], camera_filter='V', observer=settings['observer'],
                                   chunk_size=aavso_limit)
        filtered_compstars = None
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        logging.error(message)
        import traceback
        print(traceback.print_exc())
        logging.error(f"Exception during read_lightcurve for {star.path}")

    end = timer()
    logging.debug(f"Full lightcurve/phase: {end - start}")


def determine_period(df, star):
    period: Period = None
    starfile = star.get_metadata("STARFILE")
    vsx_metadata = star.get_metadata("VSX")
    is_vsx = vsx_metadata is not None
    is_vsx_with_period = is_vsx and not np.isnan(vsx_metadata.extradata['Period'])
    is_selected_with_period = starfile is not None and starfile.period is not None

    if is_selected_with_period:
        period: Period = Period(star.get_metadata("STARFILE").period, "OWN")
        logging.debug(f"Using OWN period for star {star.local_id}: {period.period}")
    elif is_vsx_with_period:
        period: Period = Period(vsx_metadata.extradata['Period'], "VSX") if vsx_metadata.extradata is \
                                                                            not None else period
        logging.debug(f"Using VSX period for star {star.local_id}: {period.period}")
    else:
        period: Period = calculate_ls_period(df.copy())
        logging.debug(f"Using LS period for star {star.local_id}: {period.period}")
    logging.debug(f"Using period: {period.period} for star {star.local_id}")
    return period


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
    CHUNK = 1  # max(1, len(star_descriptions) // nr_threads*10)
    set_seaborn_style()
    pool = mp.Pool(nr_threads, maxtasksperchild=10)
    phasedir = PurePath(resultdir, phasepart)
    chartsdir = PurePath(resultdir, chartspart)
    aavsodir = PurePath(resultdir, aavso_part)
    logging.debug(f"Using {nr_threads} threads for lightcurve, phase plotting and aavso reporting.")

    if do_phase:
        trash_and_recreate_dir(phasedir)
        trash_and_recreate_dir(PurePath(phasedir, PurePath('txt')))
    if do_light or do_light_raw:
        trash_and_recreate_dir(chartsdir)
    if do_aavso:
        trash_and_recreate_dir(aavsodir)
    comp_stars_proxy = Manager().Value('comp_stars', comp_stars)

    func = partial(read_vast_lightcurves, basedir=basedir, compstarproxy=comp_stars_proxy, do_light=do_light,
                   do_light_raw=do_light_raw, do_phase=do_phase, do_aavso=do_aavso, aavso_limit=aavsolimit,
                   phasedir=phasedir, chartsdir=chartsdir, aavsodir=aavsodir)
    with tqdm.tqdm(total=len(star_descriptions), desc=desc, unit='stars') as pbar:
        for _ in pool.imap_unordered(func, star_descriptions, chunksize=CHUNK):
            pbar.update(1)
            pass
    pool.close()
    pool.join()


# def run(star_descriptions, comparison_stars, do_charts, do_phase):


# this is a unit test
if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description='munipack automation cli')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (usually the vast dir)",
                        nargs='?', required=True)
    parser.add_argument('-s', '--starlist', help="List of stars to generate", nargs='+', required=True)
    parser.add_argument('-r', '--resultdir',
                        help="The directory where all results will be written",
                        nargs='?', required=True)
    parser.add_argument('-p', '--phase', help="Generate phase charts", action="store_true")
    parser.add_argument('-i', '--lightcurve', help="Generate lightcurve charts", action="store_true")
    parser.add_argument('-a', '--aavso', help="Generate aavso reports",
                        action="store_true")
    parser.add_argument('-k', '--checkstarfile', help="The bright and stable stars used to do ensemble photometry",
                        nargs='?', required=True)
    args = parser.parse_args()
    vastdir = utils.add_trailing_slash(args.datadir)

    # 7982
    # stars = do_calibration.get_empty_star_descriptions(args.starlist)
    stars = do_calibration.get_random_star_descriptions(len(args.starlist))
    for idx, star in enumerate(stars):
        star.local_id = int(args.starlist[idx])
        star.path = f"{vastdir}out{int(star.local_id):05}.dat"

    check_stars = main_vast.read_checkstars(args.checkstarfile)
    import ucac4
    from astropy.coordinates import SkyCoord

    real_sd = [3174, 2620]
    for idx, check_star in enumerate(check_stars):
        ucacsd = ucac4.get_ucac4_star_description(check_star)
        stars.append(StarDescription(local_id=real_sd[idx], coords=ucacsd.coords,
                                     path=f"{vastdir}out{real_sd[idx]:05}.dat"))
    comp_stars = main_vast.set_comp_stars_and_ucac4(stars, None, args.checkstarfile, vastdir)

    # run(stars, comp_stars, vastdir, args.resultdir, 'phase_candidates/', 'light_candidates/',
    run(stars[:-len(real_sd)], comp_stars, vastdir, args.resultdir, 'phase_candidates/', 'light_candidates/',
        'aavso_candidates/', do_phase=args.phase, do_light=args.light, do_aavso=args.aavso, nr_threads=1,
        desc="Phase/light/aavso of candidates")
# def run(star_descriptions, comp_stars: ComparisonStars, basedir: str, resultdir: str, phasepart: str, chartspart: str,
#         aavso_part: str, do_charts=False, do_phase=True, do_aavso=False, aavsolimit=None, nr_threads=cpu_count(),
#         desc="Writing light curve charts/phase diagrams"):

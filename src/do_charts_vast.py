import toml

import main_vast
import do_compstars
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
import tqdm
import numpy as np
import pandas as pd
import do_aavso_report
import do_calibration
import reading
import utils
import argparse
import logging
import gc
import subprocess
import math
import matplotlib as mplotlib
from multiprocessing import cpu_count, Manager
from comparison_stars import ComparisonStars
from functools import partial
from gatspy.periodic import LombScargleFast
from gatspy.periodic import TrendedLombScargle
from reading import trash_and_recreate_dir
from reading import file_selector
from timeit import default_timer as timer
from star_description import StarDescription
from pathlib import PurePath
from typing import Tuple
from pandas import DataFrame, Series

from star_metadata import StarFileData

gc.enable()
mplotlib.use('Agg')  # needs no X server
TITLE_PAD = 40


def set_seaborn_style():
    sns.set_context("notebook", font_scale=4)
    sns.set_style("ticks")


def plot_lightcurve_jd(star: StarDescription, curve: DataFrame, chartsdir, period: float):
    try:
        star_id = star.local_id
        logging.debug(f"Plotting lightcurve for {star_id}")
        vsx_name, separation, filename_no_ext = get_star_or_catalog_name(star, suffix=f"_light")
        vsx_title = '' if vsx_name is None else f" ({vsx_name} - dist:{separation:.4f})"
        save_location = PurePath(chartsdir, filename_no_ext + '.png')
        start = timer()
        upsilon_match = star.get_metadata('UPSILON') if star.has_metadata('UPSILON') else None
        upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
        end = timer()
        logging.debug(f"timing upsilon stuff {end - start}")
        coord = star.coords
        # print(f'Plotting lightcurve with star_description:{star_description}, curve length:{len(curve)}, star:{star}, curve:{curve}')
        if (curve is None):
            logging.info(f"Curve is None for star {star_id}")
            return

        curve_max = curve['realV'].max()
        curve_min = curve['realV'].min()
        curve['phaseJD'], xmin, ymin = phase_lock_lightcurve(curve['floatJD'], period)
        fig = plt.figure(figsize=(20, 12), dpi=150, facecolor='w', edgecolor='k')
        # plt.scatter(curve['floatJD'], curve['realV'])
        plt.errorbar(curve['phaseJD'], curve['realV'], yerr=curve['realErr'], linestyle='none', marker='o',
                     ecolor='gray',
                     elinewidth=1)

        plt.xlabel('phase-adjusted JD', labelpad=TITLE_PAD)
        #plt.xscale('log')
        # plt.ylabel("Absolute Mag, comp star = {:2.2f}".format(comparison_stars[0].vmag), labelpad=TITLE_PAD)
        plot_max = curve_max + 0.1
        plot_min = curve_min - 0.1
        plt.ylim(plot_min, plot_max)
        plt.xlim(xmin, ymin)
        plt.gca().invert_yaxis()
        # plt.ticklabel_format(style='plain', axis='x')
        plt.title(f"Star {star_id}{vsx_title}\nposition: {utils.get_hms_dms_matplotlib(coord)}{upsilon_text}",
                  pad=TITLE_PAD)
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


def phase_lock_lightcurve(series: Series, period):
    jds = series.to_numpy()
    jds_norm_orig = np.subtract(jds, jds.min())
    jds_norm = np.diff(jds_norm_orig)
    jds_norm = np.mod(jds_norm, period)
    jds_norm = np.concatenate(([jds_norm_orig[0]],jds_norm))
    jds_norm = np.cumsum(jds_norm)
    return jds_norm, np.min(jds_norm), np.max(jds_norm)


def plot_lightcurve(star: StarDescription, curve: DataFrame, chartsdir):
    try:
        star_id = star.local_id
        logging.debug(f"Plotting lightcurve for {star_id}")
        vsx_name, separation, filename_no_ext = get_star_or_catalog_name(star, suffix=f"_light")
        vsx_title = '' if vsx_name is None else f" ({vsx_name} - dist:{separation:.4f})"
        save_location = PurePath(chartsdir, filename_no_ext + '.png')
        start = timer()
        upsilon_match = star.get_metadata('UPSILON') if star.has_metadata('UPSILON') else None
        upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
        end = timer()
        logging.debug(f"timing upsilon stuff {end - start}")
        coord = star.coords
        # print(f'Plotting lightcurve with star_description:{star_description}, curve length:{len(curve)}, star:{star}, curve:{curve}')
        if (curve is None):
            logging.info(f"Curve is None for star {star_id}")
            return

        curve_max = curve['realV'].max()
        curve_min = curve['realV'].min()
        fig = plt.figure(figsize=(20, 12), dpi=150, facecolor='w', edgecolor='k')
        # plt.scatter(curve['floatJD'], curve['realV'])
        plt.errorbar(curve['floatJD'], curve['realV'], yerr=curve['realErr'], linestyle='none', marker='o',
                     ecolor='gray',
                     elinewidth=1)

        plt.xlabel('JD', labelpad=TITLE_PAD)
        # plt.ylabel("Absolute Mag, comp star = {:2.2f}".format(comparison_stars[0].vmag), labelpad=TITLE_PAD)
        plot_max = curve_max + 0.1
        plot_min = curve_min - 0.1
        plt.ylim(plot_min, plot_max)
        plt.xlim(curve['floatJD'].min(), curve['floatJD'].max())
        plt.gca().invert_yaxis()
        plt.ticklabel_format(style='plain', axis='x')
        plt.title(f"Star {star_id}{vsx_title}\nposition: {utils.get_hms_dms_matplotlib(coord)}{upsilon_text}",
                  pad=TITLE_PAD)
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


def plot_phase_diagram(star: StarDescription, curve: DataFrame, fullphasedir, suffix='', period: float = None):
    try:
        logging.debug(f"Starting plot phase diagram with {star} and {fullphasedir}")
        coords = star.coords
        star_id = star.local_id
        vsx_name, _, filename_no_ext = get_star_or_catalog_name(star, suffix=f"_phase{suffix}")
        vsx_title = f"{vsx_name}\n" if vsx_name is not None else ''
        save_location = PurePath(fullphasedir, filename_no_ext + '.png')
        upsilon_match = star.get_metadata('UPSILON')
        upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
        # print("Calculating phase diagram for", star)
        if curve is None:
            logging.info("Curve of star {} is None".format(star_id))
            return
        t_np = curve['floatJD']
        y_np = curve['realV'].to_numpy()
        dy_np = curve['realErr'].to_numpy()
        if period is None:
            period_max = np.max(t_np) - np.min(t_np)
            if period_max <= 0.01:
                return
            ls = LombScargleFast(optimizer_kwds={'quiet': True, 'period_range': (0.01, period_max)},
                                 silence_warnings=True).fit(t_np, y_np)
            period = ls.best_period
            # TODO test this trended lombscargle !
            # tmodel = TrendedLombScargle(optimizer_kwds={'quiet': True, 'period_range': (0.01, period_max)},
            #                             silence_warnings=True).fit(t_np, y_np)
            # period = tmodel.best_period
        fig = plt.figure(figsize=(18, 16), dpi=80, facecolor='w', edgecolor='k')
        plt.xlabel("Phase", labelpad=TITLE_PAD)
        plt.ylabel("Magnitude", labelpad=TITLE_PAD)
        plt.title(f"{vsx_title}Star {star_id}, p: {period:.5f} d{upsilon_text}\n{utils.get_hms_dms_matplotlib(coords)}",
                  pad=TITLE_PAD)
        plt.ticklabel_format(style='plain', axis='x')
        # plt.title(f"Star {star} - {period}", pad=TITLE_PAD)
        plt.tight_layout()
        # plotting + calculation of 'double' phase diagram from -1 to 1
        phased_t = np.fmod(t_np / period, 1)
        minus_one = np.vectorize(lambda t: t - 1)
        phased_t2 = minus_one(phased_t)
        phased_lc = y_np[:]
        phased_t_final = np.append(phased_t2, phased_t)
        phased_lc_final = np.append(phased_lc, phased_lc)
        phased_err = np.clip(np.append(dy_np, dy_np), -0.5, 0.5)  # error values are clipped to +0.5 and -0.5
        plt.gca().invert_yaxis()
        plt.errorbar(phased_t_final, phased_lc_final, yerr=phased_err, linestyle='none', marker='o', ecolor='gray',
                     elinewidth=1)
        logging.debug(f"Saving phase plot to {save_location}")
        fig.savefig(save_location, format='png')
        plt.close(fig)
        plt.clf()

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
            print(starfiledata, starfiledata.var_type)
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
        tomldict['period'] = float(period)
        tomldict['min'] = float(ymin)
        tomldict['max'] = float(ymax)
        tomldict['range'] = var_range
        tomldict['coords'] = [coords.ra.deg, coords.dec.deg]
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
    return period


def get_star_or_catalog_name(star: StarDescription, suffix: str):
    if star.has_metadata("VSX"):
        catalog_name, separation = star.get_metadata("VSX").get_name_and_separation()
    elif star.has_metadata("OWNCATALOG"):
        catalog_name, separation = star.get_metadata("OWNCATALOG").get_name_and_separation()
    else:
        catalog_name, separation = None, None
    filename_no_ext = f"{catalog_name}{suffix}" if catalog_name is not None else f"{star.local_id:05}{suffix}"
    return catalog_name, separation, filename_no_ext.replace(' ', '_')


def format_date(x, pos=None):
    thisind = np.clip(int(x + 0.5), 0, N - 1)
    return r.date[thisind].strftime('%Y-%m-%d')


def read_vast_lightcurves(star: StarDescription, compstarproxy, do_charts, do_phase, do_aavso,
                          aavso_limit, basedir: str, chartsdir: PurePath, phasedir: PurePath, aavsodir: PurePath):
    start = timer()
    if star.path is '':
        logging.debug(f"Path for {star.local_id} is empty")
        return
    if not do_charts and not do_phase:
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
        filtered_compstars = None
        df['floatJD'] = df['JD'].astype(np.float)
        old_size = len(df)
        # remove errors
        # df = utils.reject_outliers_iqr(df, 'realV', 20)
        # logging.info(f"Rejected {old_size-len(df)} observations with iqr.")
        if do_phase:
            if star.has_metadata("STARFILE") and star.get_metadata("STARFILE") is None:
                # This should not happen!!!
                logging.warning(f"This should not happen !!!!!!!!!!!!!!!!!!!! "
                                f"{star.local_id}, {star.get_metadata_list()}, {star.has_metadata('STARFILE')}, "
                                f"{[star.get_metadata(x) for x in star.get_metadata_list()]}")

            if star.get_metadata("STARFILE") is not None:
                period: float = star.get_metadata("STARFILE").period
                logging.debug(f"Using provided period for star {star.local_id}: {period}")
                period = plot_phase_diagram(star, df.copy(), phasedir, period=period)
            else:
                period = plot_phase_diagram(star, df.copy(), phasedir, period=None, suffix="")
        if do_charts:
            # start = timer()
            # logging.debug("NO LIGHTCRUVEGYET ")
            plot_lightcurve_jd(star, df.copy(), chartsdir, period)
            # end = timer()
            # print("plotting lightcurve:", end-start)
        if do_aavso:
            # TODO put this in settings.txt
            sitelat = '-22 57 10'
            sitelong = '-68 10 49'
            sitealt = 2500
            observer = 'ZZZ'
            do_aavso_report.report(star, df.copy(), target_dir=aavsodir, vastdir=basedir, sitelat=sitelat,
                                   sitelong=sitelong, sitealt=sitealt, filter='V',
                                   observer=observer, chunk_size=aavso_limit)
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        logging.error(message)
        import traceback
        print(traceback.print_exc())
        logging.error(f"Exception during read_lightcurve for {star.path}")

    end = timer()
    logging.debug(f"Full lightcurve/phase: {end - start}")


# reads lightcurves and passes them to lightcurve plot or phase plot
def run(star_descriptions, comp_stars: ComparisonStars, basedir: str, resultdir: str, phasepart: str, chartspart: str,
        aavso_part: str, do_charts=False, do_phase=True, do_aavso=False, aavsolimit=None, nr_threads=cpu_count(),
        desc="Writing light curve charts/phase diagrams"):
    CHUNK = 1  # max(1, len(star_descriptions) // nr_threads*10)
    set_seaborn_style()
    pool = mp.Pool(nr_threads, maxtasksperchild=10)
    phasedir = PurePath(resultdir, phasepart)
    chartsdir = PurePath(resultdir, chartspart)
    aavsodir = PurePath(resultdir, aavso_part)
    logging.info(f"Using {nr_threads} threads for lightcurve, phase plotting and aavso reporting.")

    if do_phase:
        trash_and_recreate_dir(phasedir)
        trash_and_recreate_dir(PurePath(phasedir, PurePath('txt')))
    if do_charts:
        trash_and_recreate_dir(chartsdir)
    if do_aavso:
        trash_and_recreate_dir(aavsodir)
    comp_stars_proxy = Manager().Value('comp_stars', comp_stars)

    func = partial(read_vast_lightcurves, basedir=basedir, compstarproxy=comp_stars_proxy, do_charts=do_charts,
                   do_phase=do_phase, do_aavso=do_aavso, aavso_limit=aavsolimit,
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
    comp_stars = main_vast.read_comparison_stars(stars, args.checkstarfile, vastdir)

    # run(stars, comp_stars, vastdir, args.resultdir, 'phase_candidates/', 'light_candidates/',
    run(stars[:-len(real_sd)], comp_stars, vastdir, args.resultdir, 'phase_candidates/', 'light_candidates/',
        'aavso_candidates/', do_phase=args.phase, do_charts=args.light, do_aavso=args.aavso, nr_threads=1,
        desc="Phase/light/aavso of candidates")
# def run(star_descriptions, comp_stars: ComparisonStars, basedir: str, resultdir: str, phasepart: str, chartspart: str,
#         aavso_part: str, do_charts=False, do_phase=True, do_aavso=False, aavsolimit=None, nr_threads=cpu_count(),
#         desc="Writing light curve charts/phase diagrams"):

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
import utils
import argparse
import logging
import subprocess
import math
import matplotlib as mplotlib
from multiprocessing import cpu_count
from comparison_stars import ComparisonStars
from functools import partial
from gatspy.periodic import LombScargleFast
from reading import trash_and_recreate_dir
from reading import file_selector
from timeit import default_timer as timer
from star_description import StarDescription
from pathlib import PurePath
from typing import Tuple
from pandas import DataFrame

mplotlib.use('Agg')  # needs no X server
TITLE_PAD = 40


def set_seaborn_style():
    sns.set_context("notebook", font_scale=4)
    sns.set_style("ticks")


# takes:  [ {'id': star_id, 'match': {'name': match_name, 'separation': separation_deg  } } ]
def plot_lightcurve(star_tuple: Tuple[StarDescription, DataFrame], chartsdir):
    #    try:
    star_description = star_tuple[0]
    curve = star_tuple[1]
    star_id = star_description.local_id
    logging.debug(f"Plotting light curve for {star_id}")
    vsx_name, separation, filename_no_ext = get_star_or_vsx_name(star_description, suffix=f"_light")
    vsx_title = '' if vsx_name is None else f" ({vsx_name} - dist:{separation:.4f})"
    save_location = PurePath(chartsdir, filename_no_ext)
    start = timer()
    upsilon_match = star_description.get_catalog('Upsilon')
    upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
    end = timer()
    logging.debug(f"timing upsilon stuff {end - start}")
    coord = star_description.coords
    # print(f'Plotting lightcurve with star_description:{star_description}, curve length:{len(curve)}, star:{star}, curve:{curve}')
    if (curve is None):
        logging.info(f"Curve is None for star {star_id}")
        return
    # curve = curve.replace(to_replace=99.99999, value=np.nan, inplace=False) # we filter now
    used_curve = curve
    used_curve_max = used_curve['realV'].max()
    used_curve_min = used_curve['realV'].min()
    # print(f"used curve:{used_curve['V-C']}, used curve max:{used_curve_max}, used curve min:{used_curve_min}")

    # insert counting column
    used_curve.insert(0, 'Count', range(0, len(used_curve)))
    g = sns.lmplot('Count', 'Vrel',
                   data=used_curve, height=20, aspect=5, scatter_kws={"s": 15},
                   fit_reg=False)

    plt.xlabel('Count', labelpad=TITLE_PAD)
    # plt.ylabel("Absolute Mag, comp star = {:2.2f}".format(comparison_stars[0].vmag), labelpad=TITLE_PAD)
    plot_max = used_curve_max
    plot_min = min(plot_max - 1, used_curve_min)
    logging.debug(f'min {plot_min} max {plot_max} usedmin {used_curve_min} usedmax {used_curve_max}')
    if np.isnan(plot_max) or np.isnan(plot_min):
        logging.info(f"star is nan:{star_id}, plot_max:{plot_max}, plot_min:{plot_min}")
        return
    plt.ylim(plot_min, plot_max)
    plt.xlim(0, len(used_curve))
    plt.gca().invert_yaxis()
    # g.map(plt.errorbar, 'Count', 'V-C', yerr='s1', fmt='o')
    # plt.ticklabel_format(style='plain', axis='x')
    plt.title(f"Star {star_id}{vsx_title}, position: {get_hms_dms(coord)}{upsilon_text}", pad=TITLE_PAD)
    start = timer()
    figure = g.fig
    figure.savefig(save_location)
    # g.savefig(chartsdir+str(star).zfill(5)+'_plot')
    end = timer()
    logging.debug(f"timing saving fig {end - start}")
    plt.close(g.fig)


#    except:
#        print("error", tuple)


def plot_lightcurve_jd(star_tuple: Tuple[StarDescription, DataFrame], chartsdir):
    star_description = star_tuple[0]
    curve = star_tuple[1]
    star_id = star_description.local_id
    logging.debug(f"Plotting lightcurve for {star_id}")
    vsx_name, separation, filename_no_ext = get_star_or_vsx_name(star_description, suffix=f"_light")
    vsx_title = '' if vsx_name is None else f" ({vsx_name} - dist:{separation:.4f})"
    save_location = PurePath(chartsdir, filename_no_ext+'.png')
    start = timer()
    upsilon_match = star_description.get_catalog('Upsilon')
    upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
    end = timer()
    logging.debug(f"timing upsilon stuff {end - start}")
    coord = star_description.coords
    # print(f'Plotting lightcurve with star_description:{star_description}, curve length:{len(curve)}, star:{star}, curve:{curve}')
    if (curve is None):
        logging.info(f"Curve is None for star {star_id}")
        return
    # curve = curve.replace(to_replace=99.99999, value=np.nan, inplace=False) # we filter now
    used_curve = curve
    used_curve_max = used_curve['realV'].max()
    used_curve_min = used_curve['realV'].min()

    # convert JD to int
    used_curve['realJD'] = used_curve['JD'].astype(np.float) #.astype(np.uint64)
    # fig = plt.figure(dpi=80, facecolor='w', edgecolor='k')
    fig = plt.figure(figsize=(20, 12), dpi=150, facecolor='w', edgecolor='k')
    # g = sns.lmplot('realJD', 'realV',
    #                data=used_curve, aspect=2, scatter_kws={"s": 30},
    #                fit_reg=False)
    plt.scatter(used_curve['realJD'], used_curve['realV'])
    plt.xlabel('JD', labelpad=TITLE_PAD)
    # plt.ylabel("Absolute Mag, comp star = {:2.2f}".format(comparison_stars[0].vmag), labelpad=TITLE_PAD)
    plot_max = used_curve_max
    plot_min = min(plot_max - 1, used_curve_min - 0.2)
    logging.debug(f'min {plot_min} max {plot_max} usedmin {used_curve_min} usedmax {used_curve_max}')
    if np.isnan(plot_max) or np.isnan(plot_min):
        logging.info(f"star is nan:{star_id}, plot_max:{plot_max}, plot_min:{plot_min}")
        return
    plt.ylim(plot_min, plot_max)
    plt.xlim(used_curve['realJD'].min(), used_curve['realJD'].max())
    plt.gca().invert_yaxis()
    # g.map(plt.errorbar, 'Count', 'V-C', yerr='s1', fmt='o')
    plt.ticklabel_format(style='plain', axis='x')
    plt.title(f"Star {star_id}{vsx_title}\nposition: {get_hms_dms(coord)}{upsilon_text}", pad=TITLE_PAD)
    plt.tight_layout()
    # figure = g.fig
    # figure.suptitle(f"Star {star_id}{vsx_title}, position: {get_hms_dms(coord)}{upsilon_text}")
    start = timer()

    fig.savefig(save_location)
    # g.savefig(chartsdir+str(star).zfill(5)+'_plot')
    end = timer()
    logging.debug(f"timing saving fig {end - start}")
    plt.close(fig)
    plt.clf()


def plot_phase_diagram(star_tuple: Tuple[StarDescription, DataFrame], fullphasedir, suffix='', period=None):
    logging.debug(f"Starting plot phase diagram with {star_tuple} and {fullphasedir}")
    star_description = star_tuple[0]
    coords = star_description.coords
    curve = star_tuple[1]
    star_id = star_description.local_id
    vsx_name, _, filename_no_ext = get_star_or_vsx_name(star_description, suffix=f"_phase{suffix}")
    vsx_title = f"{vsx_name}\n" if vsx_name is not None else ''
    save_location = PurePath(fullphasedir, filename_no_ext+'.png')
    upsilon_match = star_description.get_catalog('Upsilon')
    upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
    # print("Calculating phase diagram for", star)
    if curve is None:
        logging.info("Curve of star {} is None".format(star_id))
        return
    t_np = curve['JD'].astype(np.float)
    y_np = curve['realV'].to_numpy()
    dy_np = curve['realErr'].to_numpy()
    if period is None:
        period_max = np.max(t_np) - np.min(t_np)
        if period_max <= 0.01:
            return
        ls = LombScargleFast(optimizer_kwds={'quiet': True, 'period_range': (0.01, period_max)}, silence_warnings=True) \
            .fit(t_np, y_np)
        period = ls.best_period
    fig = plt.figure(figsize=(18, 16), dpi=80, facecolor='w', edgecolor='k')
    plt.xlabel("Phase", labelpad=TITLE_PAD)
    plt.ylabel("Magnitude", labelpad=TITLE_PAD)
    plt.title(f"{vsx_title}Star {star_id}, p: {period:.5f} d{upsilon_text}\n{get_hms_dms(coords)}", pad=TITLE_PAD)
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
    with open(f"{fullphasedir}/txt/{filename_no_ext}.txt", 'w') as f:
        f.write('\n'.join([f"period={period}", f'range="{np.min(y_np):.1f}-{np.max(y_np):.1f}"',
                           f"coords=[{coords.ra.deg}, {coords.dec.deg}]"]))


def get_star_or_vsx_name(star_description: StarDescription, suffix: str):
    vsx_name, separation = star_description.get_match_string("VSX")
    filename_no_ext = f"{vsx_name}{suffix}" if vsx_name is not None else f"{star_description.local_id:05}{suffix}"
    return vsx_name, separation, filename_no_ext.replace(' ', '_')


def reject_outliers_iqr(data, cut=5):
    q1, q3 = np.percentile(data, [cut, 100 - cut])
    iqr = q3 - q1
    lower_bound = q1 - (iqr * 1.5)
    upper_bound = q3 + (iqr * 1.5)
    return np.where((data > lower_bound) & (data < upper_bound))


def get_hms_dms(coord):
    return "{:2.0f}$^h$ {:02.0f}$^m$ {:02.2f}$^s$ | {:2.0f}$\degree$ {:02.0f}$'$ {:02.2f}$''$" \
        .format(coord.ra.hms.h, abs(coord.ra.hms.m), abs(coord.ra.hms.s),
                coord.dec.dms.d, abs(coord.dec.dms.m), abs(coord.dec.dms.s))


def format_date(x, pos=None):
    thisind = np.clip(int(x + 0.5), 0, N - 1)
    return r.date[thisind].strftime('%Y-%m-%d')


def read_vast_lightcurves(star_description: StarDescription, comp_stars: ComparisonStars, do_charts, do_phase, do_aavso,
                          aavso_limit, basedir: str, chartsdir: PurePath, phasedir: PurePath, aavsodir: PurePath):
    start = timer()
    if star_description.path is '':
        logging.debug(f"Path for {star_description.local_id} is empty")
        return
    if not do_charts and not do_phase:
        logging.debug("Nothing to do, no charts or phase needed")

    logging.debug(
        f"Reading lightcurves for star {star_description.local_id} at path {star_description.path} for {star_description}...")
    # comp_mags = [x.vmag for x in comparison_stars]

    try:
        df = pd.read_csv(star_description.path, delim_whitespace=True,
                         names=['JD', 'Vrel', 'err', 'X', 'Y', 'unknown', 'file', 'vast1', 'vast2', 'vast3', 'vast4',
                                'vast5', 'vast6', 'vast7', 'vast8', 'vast9', 'vast10'], dtype={'JD': str})

        if df is None or len(df) == 0:
            logging.info(f"No lightcurve found for {star_description.path}")
            return

        # adding vmag of comparison star to all diff mags
        # TODO: this is wrong, should be composition of all comparison stars. To do error propagation use quadrature
        # rule: http://ipl.physics.harvard.edu/wp-uploads/2013/03/PS3_Error_Propagation_sp13.pdf
        # df['V-C'] = df['V-C'] + comparison_stars[0].vmag
        close_comp_stars = do_compstars.closest_compstars(star_description, comp_stars, 10)
        df['realV'], df['realErr'] = do_compstars.calculate_mean_value_ensemble_photometry(df, close_comp_stars)
        tuple = star_description, df
        if do_charts:
            # start = timer()
            # logging.debug("NO LICGHTCRUVEGYET ")
            plot_lightcurve_jd(tuple, chartsdir)
            # end = timer()
            # print("plotting lightcurve:", end-start)
        if do_phase:
            # start = timer()
            try:
                phase = subprocess.check_output([f"{basedir}lib/BLS/bls", star_description.path])
            except:
                phase = None

            if phase is not None:
                plot_phase_diagram(tuple, phasedir, period=float(phase.decode("utf-8")[:-2]))
            plot_phase_diagram(tuple, phasedir, period=None, suffix="")
            # end = timer()
            # print("plotting phase:", end-start)
        if do_aavso:
            # TODO put this in settings.txt
            sitelat = '-22 57 10'
            sitelong = '-68 10 49'
            sitealt = 2500
            observer = 'ZZZ'

            do_aavso_report.report(tuple, target_dir=aavsodir, vastdir=basedir, sitelat=sitelat, sitelong=sitelong,
                                   sitealt=sitealt, comparison_stars=comp_stars, filter='V',
                                   observer=observer, chunk_size=aavso_limit)

    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        logging.error(message)
        logging.error("File not found error in store and curve for star", star_description.path)

    end = timer()
    logging.debug(f"Full lightcurve/phase: {end - start}")


# reads lightcurves and passes them to lightcurve plot or phase plot
def run(star_descriptions, comp_stars: ComparisonStars, basedir: str, resultdir: str, phasepart: str, chartspart: str,
        aavso_part: str, do_charts=False, do_phase=True, do_aavso=False, aavsolimit=None, nr_threads=cpu_count(),
        desc="Writing light curve charts/phase diagrams"):
    CHUNK = 1
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

    func = partial(read_vast_lightcurves, basedir=basedir, comp_stars=comp_stars, do_charts=do_charts,
                   do_phase=do_phase, do_aavso=do_aavso, aavso_limit=aavsolimit,
                   phasedir=phasedir, chartsdir=chartsdir, aavsodir=aavsodir)
    with tqdm.tqdm(total=len(star_descriptions), desc=desc) as pbar:
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

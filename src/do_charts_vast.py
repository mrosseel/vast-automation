from functools import partial
from multiprocessing import cpu_count
from comparison_stars import ComparisonStars
import reading
import matplotlib as mp

mp.use('Agg')  # needs no X server
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
import tqdm
import numpy as np
import pandas as pd
from gatspy.periodic import LombScargleFast
from reading import trash_and_recreate_dir
from reading import file_selector
import do_calibration
from timeit import default_timer as timer
import utils
import argparse
from star_description import StarDescription
import logging
import subprocess
import math

TITLE_PAD = 40


def set_seaborn_style():
    sns.set_context("notebook", font_scale=4)
    sns.set_style("ticks")


# takes:  [ {'id': star_id, 'match': {'name': match_name, 'separation': separation_deg  } } ]
def plot_lightcurve(tuple, chartsdir):
    #    try:
    star_description = tuple[0]
    curve = tuple[1]
    star = star_description.local_id
    logging.debug(f"Plotting lightcurve for {star}")
    star_match, separation = star_description.get_match_string("VSX")
    match_string = f"({star_match})" if not star_match == None else ''
    star_name = '' if star_match == None else " ({} - dist:{:.4f})".format(star_match, separation)
    start = timer()
    upsilon_match = star_description.get_catalog('Upsilon')
    upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
    end = timer()
    logging.debug(f"timing upsilon stuff {end - start}")
    coord = star_description.coords
    # print(f'Plotting lightcurve with star_description:{star_description}, curve length:{len(curve)}, star:{star}, curve:{curve}')
    if (curve is None):
        logging.info(f"Curve is None for star {star}")
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
        logging.info(f"star is nan:{star}, plot_max:{plot_max}, plot_min:{plot_min}")
        return
    plt.ylim(plot_min, plot_max)
    plt.xlim(0, len(used_curve))
    plt.gca().invert_yaxis()
    # g.map(plt.errorbar, 'Count', 'V-C', yerr='s1', fmt='o')
    # plt.ticklabel_format(style='plain', axis='x')
    plt.title("Star {0}{1}, position: {2}{3}".format(star, star_name, get_hms_dms(coord), upsilon_text), pad=TITLE_PAD)
    start = timer()
    figure = g.fig
    figure.savefig(chartsdir + str(star).zfill(5) + '_plot')
    # g.savefig(chartsdir+str(star).zfill(5)+'_plot')
    end = timer()
    logging.debug(f"timing saving fig {end - start}")
    plt.close(g.fig)


#    except:
#        print("error", tuple)


def plot_lightcurve_jd(tuple, chartsdir):
    star_description = tuple[0]
    curve = tuple[1]
    star = star_description.local_id
    logging.debug(f"Plotting lightcurve for {star}")
    star_match, separation = star_description.get_match_string("VSX")
    match_string = f"({star_match})" if not star_match == None else ''
    star_name = '' if star_match == None else " ({} - dist:{:.4f})".format(star_match, separation)
    start = timer()
    upsilon_match = star_description.get_catalog('Upsilon')
    upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
    end = timer()
    logging.debug(f"timing upsilon stuff {end - start}")
    coord = star_description.coords
    # print(f'Plotting lightcurve with star_description:{star_description}, curve length:{len(curve)}, star:{star}, curve:{curve}')
    if (curve is None):
        logging.info(f"Curve is None for star {star}")
        return
    # curve = curve.replace(to_replace=99.99999, value=np.nan, inplace=False) # we filter now
    used_curve = curve
    used_curve_max = used_curve['realV'].max()
    used_curve_min = used_curve['realV'].min()
    # print(f"used curve:{used_curve['V-C']}, used curve max:{used_curve_max}, used curve min:{used_curve_min}")

    # convert JD to float (descructive !!!)
    used_curve['JD'] = used_curve['JD'].astype(np.float)
    g = sns.lmplot('JD', 'realV',
                   data=used_curve, height=20, aspect=5, scatter_kws={"s": 15},
                   fit_reg=False)

    plt.xlabel('JD', labelpad=TITLE_PAD)
    # plt.ylabel("Absolute Mag, comp star = {:2.2f}".format(comparison_stars[0].vmag), labelpad=TITLE_PAD)
    plot_max = used_curve_max
    plot_min = min(plot_max - 1, used_curve_min)
    logging.debug(f'min {plot_min} max {plot_max} usedmin {used_curve_min} usedmax {used_curve_max}')
    if np.isnan(plot_max) or np.isnan(plot_min):
        logging.info(f"star is nan:{star}, plot_max:{plot_max}, plot_min:{plot_min}")
        return
    plt.ylim(plot_min, plot_max)
    plt.xlim(0, len(used_curve))
    plt.gca().invert_yaxis()
    # g.map(plt.errorbar, 'Count', 'V-C', yerr='s1', fmt='o')
    # plt.ticklabel_format(style='plain', axis='x')
    plt.title("Star {0}{1}, position: {2}{3}".format(star, star_name, get_hms_dms(coord), upsilon_text), pad=TITLE_PAD)
    start = timer()
    figure = g.fig
    figure.savefig(chartsdir + str(star).zfill(5) + '_plot')
    # g.savefig(chartsdir+str(star).zfill(5)+'_plot')
    end = timer()
    logging.debug(f"timing saving fig {end - start}")
    plt.close(g.fig)


def plot_phase_diagram(tuple, fullphasedir, suffix='', period=None):
    star_description = tuple[0]
    coords = star_description.coords
    curve = tuple[1]
    star = star_description.local_id
    star_match, separation = star_description.get_match_string("VSX")
    match_string = f"{star_match}\n" if not star_match == None else ''
    upsilon_match = star_description.get_catalog('Upsilon')
    upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
    # print("Calculating phase diagram for", star)
    if curve is None:
        logging.info("Curve of star {} is None".format(star))
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
    # print("Best period: " + str(period) + " days")
    fig = plt.figure(figsize=(18, 16), dpi=80, facecolor='w', edgecolor='k')
    plt.xlabel("Phase", labelpad=TITLE_PAD)
    plt.ylabel("Magnitude", labelpad=TITLE_PAD)
    plt.title(f"{match_string}Star {star}, p: {period:.5f} d{upsilon_text}\n{get_hms_dms(coords)}", pad=TITLE_PAD)
    # plt.title(f"Star {star} - {period}", pad=TITLE_PAD)
    plt.tight_layout()
    # plotting + calculation of 'double' phase diagram from -1 to 1
    phased_t = np.fmod(t_np / period, 1)
    minus_one = lambda t: t - 1
    minus_oner = np.vectorize(minus_one)
    phased_t2 = minus_oner(phased_t)
    phased_lc = y_np[:]
    phased_t_final = np.append(phased_t2, phased_t)
    phased_lc_final = np.append(phased_lc, phased_lc)
    phased_err = np.clip(np.append(dy_np, dy_np), -0.5, 0.5)  # error values are clipped to +0.5 and -0.5
    plt.gca().invert_yaxis()
    plt.errorbar(phased_t_final, phased_lc_final, yerr=phased_err, linestyle='none', marker='o', ecolor='gray',
                 elinewidth=1)
    # save_location = phasedir+str(star).zfill(5)+'_phase'+suffix
    filename_no_ext = f"{star_match}_phase{suffix}" if star_match is not None else f"{star:05}_phase{suffix}"
    save_location = f"{fullphasedir}{filename_no_ext}" if star_match is not None else f"{fullphasedir}{filename_no_ext}"
    logging.debug(f"Saving phase plot to {save_location}")
    fig.savefig(f"{save_location}.png")
    plt.close(fig)
    with open(f"{fullphasedir}/txt/{filename_no_ext}.txt", 'w') as f:
        f.write('\n'.join([f"period={period}", f'range="{np.min(y_np):.1f}-{np.max(y_np):.1f}"',
                           f"coords=[{coords.ra.deg}, {coords.dec.deg}]"]))


def get_hms_dms(coord):
    return "{:2.0f}$^h$ {:02.0f}$^m$ {:02.2f}$^s$ | {:2.0f}$\degree$ {:02.0f}$'$ {:02.2f}$''$" \
        .format(coord.ra.hms.h, abs(coord.ra.hms.m), abs(coord.ra.hms.s),
                coord.dec.dms.d, abs(coord.dec.dms.m), abs(coord.dec.dms.s))


def format_date(x, pos=None):
    thisind = np.clip(int(x + 0.5), 0, N - 1)
    return r.date[thisind].strftime('%Y-%m-%d')


def read_vast_lightcurves(star_description: StarDescription, comp_stars: ComparisonStars, do_charts, do_phase,
                          basedir: str, chartsdir: str, phasedir: str):
    start = timer()
    if star_description.path is '':
        logging.debug(f"Path for {star_description.local_id} is empty")
        return
    logging.debug(
        f"Reading lightcurves for star {star_description.local_id} at path {star_description.path} for {star_description}...")
    # comp_mags = [x.vmag for x in comparison_stars]

    try:
        df = pd.read_csv(star_description.path, delim_whitespace=True,
                         names=['JD', 'Vrel', 'err', 'X', 'Y', 'file', 'vast1', 'vast2', 'vast3', 'vast4', 'vast5',
                                'vast6', 'vast7', 'vast8', 'vast9', 'vast10', 'vast11'], dtype={'JD': str})

        if df is None or len(df) == 0:
            logging.info(f"No lightcurve found for {star_description.path}")
            return

        # adding vmag of comparison star to all diff mags
        # TODO: this is wrong, should be composition of all comparison stars. To do error propagation use quadrature
        # rule: http://ipl.physics.harvard.edu/wp-uploads/2013/03/PS3_Error_Propagation_sp13.pdf
        # df['V-C'] = df['V-C'] + comparison_stars[0].vmag
        df['realV'], df['realErr'] = calculate_real_mag_and_err(df, comp_stars)
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
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        logging.error(message)
        logging.error("File not found error in store and curve for star", star_description.path)

    end = timer()
    logging.debug(f"Full lightcurve/phase: {end - start}")


def calculate_real_mag_and_err(df, comp_stars: ComparisonStars):
    logging.debug(f"Start calculate_real with {df.shape[0]} rows and {len(comp_stars.observations)} comp stars.")

    # the average of the real magnitude of the comparison stars
    meanreal = np.mean(comp_stars.comp_catalogmags)
    logging.debug(f"meanreal is {meanreal}")
    # error = sqrt((vsig**2+(1/n sum(sigi)**2)))
    realV = []
    realErr = []
    for index, row in df.iterrows():
        logging.debug(f"row is {row}, comp_mags is {comp_stars.comp_catalogmags}")
        comp_obs = []
        comp_err = []
        for compstar in comp_stars.observations:
            jd = row['JD']
            if jd in compstar:
                comp_obs.append(float(compstar[jd][0]))
                comp_err.append(float(compstar[jd][1]))
            else:
                logging.info(f"Key error for {row['JD']}, {comp_stars.ids}")

        meanobs = -1
        meanerr = -1
        if len(comp_obs) > 0 and len(comp_obs) == len(comp_err):
            meanobs = np.mean(comp_obs)
            # Vobs - Cobs + Creal = V
            realV.append(row['Vrel']-meanobs+meanreal)
            meanerr = math.sqrt(math.pow(row['err'], 2)+math.pow(np.mean(comp_err), 2))
            realErr.append(meanerr)
        else:  # error in the comparison stars
            realV.append(row['Vrel'])
            realErr.append(row['err'])
        if meanobs == -1 or meanerr == -1:
            logging.info(f"vrel: {row['Vrel']}, meanobs: {meanobs}, compobs: {comp_obs},  meanreal: {meanreal}, "
                         f"real {comp_stars.comp_catalogmags},  vrel: {row['Vrel'] - meanobs + meanreal}, meanerr: {meanerr},"
                         f"nr of compstar observations={len(compstar)}, nr of variable observations={len(df)}")
            logging.info(f"{len(comp_obs)}, {len(comp_obs)} == {len(comp_err)}")
        logging.debug(f"Results of photometry: V diff: {df['Vrel'].mean()-np.mean(realV)}, err diff: {df['err'].mean()-np.mean(realErr)}")
    return realV, realErr


# reads lightcurves and passes them to lightcurve plot or phase plot
# def run(star_descriptions, comparison_stars, do_charts, do_phase):
def run(star_descriptions, comp_stars: ComparisonStars, basedir: str, phasedir: str, chartsdir: str,
        do_charts=False, do_phase=True, nr_threads=cpu_count()):
    CHUNK = 1
    set_seaborn_style()
    pool = mp.Pool(nr_threads)
    # pool = mp.Pool(1)
    logging.info(f"Using {nr_threads} threads for lightcurve and phase plotting.")

    if do_phase:
        trash_and_recreate_dir(phasedir)
        trash_and_recreate_dir(phasedir + '/txt')
    if do_charts:
        trash_and_recreate_dir(chartsdir)

    func = partial(read_vast_lightcurves, basedir=basedir, comp_stars=comp_stars, do_charts=do_charts, do_phase=do_phase,
                   phasedir=phasedir, chartsdir=chartsdir)
    with tqdm.tqdm(total=len(star_descriptions), desc='Writing light curve charts/phase diagrams') as pbar:
        for _ in pool.imap_unordered(func, star_descriptions, chunksize=CHUNK):
            pbar.update(1)
            pass


# this is a unit test
if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description='munipack automation cli')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (usually the vast dir)",
                        nargs='?', required=True)
    parser.add_argument('-s', '--starlist', help="List of stars to generate", nargs='+')
    args = parser.parse_args()
    vastdir = utils.add_trailing_slash(args.datadir)

    # 7982
    stars = do_calibration.get_empty_star_descriptions(args.starlist)
    for star in stars:
        star.path = f"{vastdir}out{int(star.local_id):05}.dat"
    run(stars, "compstars", vastdir, 'phase_starlist/', 'chart_starlist/', nr_threads=cpu_count())

from functools import partial

from init_loader import init, settings
import reading
import matplotlib as mp
mp.use('Agg') # needs no X server
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
import tqdm
import numpy as np
import pandas as pd
from gatspy.periodic import LombScargleFast
from init_loader import init, settings
from reading import trash_and_recreate_dir
import logging
from timeit import default_timer as timer

TITLE_PAD=40

def set_seaborn_style():
    sns.set_context("notebook", font_scale=4)
    sns.set_style("ticks")

# takes:  [ {'id': star_id, 'match': {'name': match_name, 'separation': separation_deg  } } ]
def plot_lightcurve(tuple, comparison_stars):
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
    logging.debug(f"timing upsilon stuff {end-start}")
    coord = star_description.coords
    #print(f'Plotting lightcurve with star_description:{star_description}, curve length:{len(curve)}, star:{star}, curve:{curve}')
    if(curve is None):
        logging.info(f"Curve is None for star {star}")
        return
    #curve = curve.replace(to_replace=99.99999, value=np.nan, inplace=False) # we filter now
    used_curve = curve
    used_curve_max = used_curve['V-C'].max()
    used_curve_min = used_curve['V-C'].min()
    #print(f"used curve:{used_curve['V-C']}, used curve max:{used_curve_max}, used curve min:{used_curve_min}")

    #insert counting column
    used_curve.insert(0, 'Count', range(0, len(used_curve)))
    g = sns.lmplot('Count', 'V-C',
               data=used_curve, height=20, aspect=5,scatter_kws={"s": 15},
               fit_reg=False)

    plt.xlabel('Count', labelpad=TITLE_PAD)
    plt.ylabel("Absolute Mag, comp star = {:2.2f}".format(comparison_stars[0].vmag), labelpad=TITLE_PAD)
    plot_max = used_curve_max
    plot_min = min(plot_max-1, used_curve_min)
    logging.debug(f'min {plot_min} max {plot_max} usedmin {used_curve_min} usedmax {used_curve_max}')
    if np.isnan(plot_max) or np.isnan(plot_min):
        logging.info(f"star is nan:{star}, plot_max:{plot_max}, plot_min:{plot_min}")
        return
    plt.ylim(plot_min,plot_max)
    plt.xlim(0, len(used_curve))
    plt.gca().invert_yaxis()
    #g.map(plt.errorbar, 'Count', 'V-C', yerr='s1', fmt='o')
    #plt.ticklabel_format(style='plain', axis='x')
    plt.title("Star {0}{1}, position: {2}{3}".format(star, star_name, get_hms_dms(coord), upsilon_text), pad=TITLE_PAD)
    start = timer()
    figure = g.fig
    figure.savefig(settings.chartsdir+str(star).zfill(5)+'_plot')
    # g.savefig(settings.chartsdir+str(star).zfill(5)+'_plot')
    end = timer()
    logging.debug(f"timing saving fig {end-start}")
    plt.close(g.fig)
#    except:
#        print("error", tuple)


def plot_phase_diagram(tuple, suffix='', period=None):
    star_description = tuple[0]
    coords = star_description.coords
    curve = tuple[1]
    star = star_description.local_id
    star_match, separation = star_description.get_match_string("VSX")
    match_string = f" ({star_match})" if not star_match == None else ''
    upsilon_match = star_description.get_catalog('Upsilon')
    upsilon_text = upsilon_match.get_upsilon_string() if upsilon_match is not None else ''
    #print("Calculating phase diagram for", star)
    if curve is None:
        logging.info("Curve of star {} is None".format(star))
        return
    t_np = curve['JD'].to_numpy()
    y_np = curve['realV'].to_numpy()
    dy_np = curve['s1'].to_numpy()
    if period is None:
        period_max = np.max(t_np)-np.min(t_np)
        if period_max <= 0.01:
             return
        ls = LombScargleFast(optimizer_kwds={'quiet': True, 'period_range': (0.01,period_max)}, silence_warnings=True)\
            .fit(t_np,y_np)
        period = ls.best_period
    #print("Best period: " + str(period) + " days")
    fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    plt.xlabel("Phase", labelpad=TITLE_PAD)
    plt.ylabel("Magnitude", labelpad=TITLE_PAD)
    plt.title(f"Star {star}{match_string}, p: {period:.5f} d{upsilon_text}\n{get_hms_dms(coords)}", pad=TITLE_PAD)
    plt.tight_layout()
    # plotting + calculation of 'double' phase diagram from -1 to 1
    phased_t = np.fmod(t_np/period,1)
    minus_one = lambda t: t - 1
    minus_oner = np.vectorize(minus_one)
    phased_t2 = minus_oner(phased_t)
    phased_lc = y_np[:]
    phased_t_final = np.append(phased_t2, phased_t)
    phased_lc_final = np.append(phased_lc, phased_lc)
    phased_err = np.clip(np.append(dy_np, dy_np), -0.5, 0.5) # error values are clipped to +0.5 and -0.5
    plt.gca().invert_yaxis()
    plt.errorbar(phased_t_final,phased_lc_final,yerr=phased_err,linestyle='none',marker='o', ecolor='gray', elinewidth=1)
    fig.savefig(settings.phasedir+str(star).zfill(5)+'_phase'+suffix)
    plt.close(fig)

def get_hms_dms(coord):
    return "{:2.0f}$^h$ {:02.0f}$^m$ {:02.2f}$^s$ | {:2.0f}$\degree$ {:02.0f}$'$ {:02.2f}$''$"\
        .format(coord.ra.hms.h, abs(coord.ra.hms.m), abs(coord.ra.hms.s),
                coord.dec.dms.d, abs(coord.dec.dms.m), abs(coord.dec.dms.s))

def format_date(x, pos=None):
    thisind = np.clip(int(x + 0.5), 0, N - 1)
    return r.date[thisind].strftime('%Y-%m-%d')

def read_lightcurves(star_description, comparison_stars, do_charts, do_phase):
    start = timer()
    logging.debug("Reading lightcurves...")
    comp_mags = [x.vmag for x in comparison_stars]

    try:
        df = reading.read_lightcurve(star_description.local_id,filter=True)
        if df is None or len(df) == 0:
            logging.info(f"No lightcurve found for star {star_description.local_id}")
            return

        # adding vmag of comparison star to all diff mags
        # TODO: this is wrong, should be composition of all comparison stars. To do error propagation use quadrature
        # rule: http://ipl.physics.harvard.edu/wp-uploads/2013/03/PS3_Error_Propagation_sp13.pdf
        #df['V-C'] = df['V-C'] + comparison_stars[0].vmag
        logging.debug(f"Calculate real v for {star_description.local_id}")
        df['realV'] = calculate_real_mag(df, comp_mags)
        logging.debug(f"charting {star_description}")
        tuple = star_description, df
        if do_charts:
            # start = timer()
            plot_lightcurve(tuple, comparison_stars=comparison_stars)
            # end = timer()
            # print("plotting lightcurve:", end-start)
        if do_phase:
            # start = timer()
            plot_phase_diagram(tuple)
            # end = timer()
            # print("plotting phase:", end-start)
    except FileNotFoundError:
        logging.error("File not found error in store and curve for star", star_description.local_id)

    end = timer()
    logging.debug(f"Full lightcurve/phase: {end-start}")

def calculate_real_mag(df, comp_mags):
    logging.debug(f"Start calculate_real with {df.shape[0]} rows")

    cees = []
    for index, row in df.iterrows():
        logging.debug(f"row is {row}, comp_mags is {comp_mags}")
        the_mask = [not bool(int(x)) for x in row['mask']] # read and invert for masking
        logging.debug(f"The mask is {the_mask}")
        comp_mags_masked = np.ma.array(comp_mags, mask=the_mask)
        logging.debug(f"The masked comps is {comp_mags_masked}")
        mean = np.ma.mean(comp_mags_masked)
        logging.debug(f"The mean is {mean}")
        cees.append(mean)
        logging.debug(f"print cees: {cees}")

    logging.debug(f"print cees: {cees}")
    logging.debug(f"df V-C beforel: {df['V-C']}")
    result  = df['V-C'].add(cees)
    logging.debug(f"df V-C after: {result}")
    logging.debug(f"end calculate_real")
    return result



# reads lightcurves and passes them to lightcurve plot or phase plot
def run(star_descriptions, comparison_stars, do_charts, do_phase):
    CHUNK = 2
    set_seaborn_style()
    pool = mp.Pool(init.nr_threads)
    logging.info(f"Using {init.nr_threads} threads for lightcurve and phase plotting.")

    if do_charts:
        trash_and_recreate_dir(settings.chartsdir)
    if do_phase:
        trash_and_recreate_dir(settings.phasedir)

    func = partial(read_lightcurves, comparison_stars=comparison_stars, do_charts=do_charts, do_phase=do_phase)
    with tqdm.tqdm(total=len(star_descriptions), desc='Writing light curve charts/phase diagrams') as pbar:
        for _ in pool.imap_unordered(func, star_descriptions, chunksize=CHUNK):
            pbar.update(1)
            pass

# this is a unit test
if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    a = np.zeros(shape=(3,2))
    df = pd.DataFrame(a,columns=['V-C','mask'])
    df['V-C'] = pd.to_numeric(df['V-C'])
    df['mask'] = ['01', '11', '10']
    result = calculate_real_mag(df, [10, 12])
    print('\n\n\n',result)
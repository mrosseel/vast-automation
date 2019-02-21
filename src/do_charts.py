from functools import partial

import init
import reading
#import astropy_helper
import matplotlib as mp
mp.use('Agg') # needs no X server
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
import tqdm
import numpy as np
from astropy.coordinates import SkyCoord
from gatspy.periodic import LombScargleFast
import init
from star_description import get_upsilon_string
from reading import trash_and_recreate_dir

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
    star_match, separation = star_description.get_match_string("VSX", strict=False)
    match_string = f"({star_match})" if not star_match == None else ''
    star_name = '' if star_match == None else " ({} - dist:{:.4f})".format(star_match, separation)
    upsilon_text = get_upsilon_string(star_description)
    coord = star_description.coords
    #print(f'Plotting lightcurve with star_description:{star_description}, curve length:{len(curve)}, star:{star}, curve:{curve}')
    if(curve is None):
        print("Curve is None for star", star)
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

    plt.title("Star {0}{1}, position: {2}{3}".format(star, star_name, get_hms_dms(coord), upsilon_text), pad=TITLE_PAD)

    #plt.xaxis.set_major_formatter(ticker.FuncFormatter(format_date))
    #plt.set_title("Custom tick formatter")
    #fig.autofmt_xdate()
    plt.xlabel('Count', labelpad=TITLE_PAD)
    plt.ylabel("Absolute Mag, comp star = {:2.2f}".format(comparison_stars[0].vmag), labelpad=TITLE_PAD)
    plot_max = used_curve_max
    plot_min = min(plot_max-1, used_curve_min)
    #print('min', plot_min, 'max', plot_max, 'usedmin', used_curve_min, 'usedmax', used_curve_max)
    if np.isnan(plot_max) or np.isnan(plot_min):
        print(f"star is nan:{star}, plot_max:{plot_max}, plot_min:{plot_min}")
        return
    #print("Star:{},dim:{},bright:{}".format(star, plot_dim, plot_bright))
    plt.ylim(plot_min,plot_max)
    plt.xlim(0, len(used_curve))
    plt.gca().invert_yaxis()
    #g.map(plt.errorbar, 'Count', 'V-C', yerr='s1', fmt='o')
    #plt.ticklabel_format(style='plain', axis='x')
    g.savefig(init.chartsdir+str(star).zfill(5)+'_plot' )
    plt.close(g.fig)
#    except:
#        print("error", tuple)

def plot_phase_diagram(tuple, comparison_stars, suffix='', period=None):
    star_description = tuple[0]
    curve = tuple[1]
    star = star_description.local_id
    star_match, separation = star_description.get_match_string("VSX", strict=False)
    match_string = f" ({star_match})" if not star_match == None else ''
    upsilon_text = get_upsilon_string(star_description)
    #print("Calculating phase diagram for", star)
    if curve is None:
        print("Curve of star {} is None".format(star))
        return
    t_np = curve['JD'].values()
    y_np = curve['V-C'].values()
    dy_np = curve['s1'].values()
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
    plt.title("Star {}{}, period: {:.5f} d{}".format(star, match_string, period, upsilon_text), pad=TITLE_PAD)
    plt.tight_layout()
    # plotting + calculation of 'double' phase diagram from -1 to 1
    phased_t = np.fmod(t_np/period,1)
    minus_one = lambda t: t - 1
    minus_oner = np.vectorize(minus_one)
    phased_t2 = minus_oner(phased_t)
    phased_lc = y_np[:]
    phased_t_final = np.append(phased_t2, phased_t)
    phased_lc_final = np.append(phased_lc, phased_lc)
    #phased_lc_final = phased_lc_final + comparison_stars[0].vmag
    phased_err = np.append(dy_np, dy_np)
    plt.gca().invert_yaxis()
    plt.errorbar(phased_t_final,phased_lc_final,yerr=phased_err,linestyle='none',marker='o', ecolor='gray', elinewidth=1)
    fig.savefig(init.phasedir+str(star).zfill(5)+'_phase'+suffix)
    plt.close(fig)

def get_hms_dms(coord):
    return "{:2.0f}$^h$ {:02.0f}$^m$ {:02.2f}$^s$ | {:2.0f}$\degree$ {:02.0f}$'$ {:02.2f}$''$"\
        .format(coord.ra.hms.h, abs(coord.ra.hms.m), abs(coord.ra.hms.s),
                coord.dec.dms.d, abs(coord.dec.dms.m), abs(coord.dec.dms.s))

def format_date(x, pos=None):
    thisind = np.clip(int(x + 0.5), 0, N - 1)
    return r.date[thisind].strftime('%Y-%m-%d')

def read_lightcurves(star_pos, star_descriptions, comparison_stars, do_charts, do_phase):
    star_description = star_descriptions[star_pos]
    try:
        df = reading.read_lightcurve(star_description.local_id,filter=True)
        if df is None or len(df) == 0:
            print("No lightcurve found for star", star_description.local_id)
            return
        # adding vmag of comparison star to all diff mags
        df['V-C'] = df['V-C'] + comparison_stars[0].vmag
        tuple = star_description, df
        if do_charts:
            plot_lightcurve(tuple, comparison_stars=comparison_stars)
        if do_phase:
            plot_phase_diagram(tuple, comparison_stars=comparison_stars)
    except FileNotFoundError:
        print("File not found error in store and curve for star", star_description.local_id)

# reads lightcurves and passes them to lightcurve plot or phase plot
def run(star_descriptions, comparison_stars, do_charts, do_phase):
    CHUNK = 100
    # star_list = [star.local_id for star in star_descriptions]
    star_list = range(0, len(star_descriptions))
    set_seaborn_style()
    pool = mp.Pool(init.nr_threads)

    if do_charts:
        trash_and_recreate_dir(init.chartsdir)
    if do_phase:
        trash_and_recreate_dir(init.phasedir)

    func = partial(read_lightcurves, star_descriptions=star_descriptions, comparison_stars=comparison_stars, do_charts=do_charts, do_phase=do_phase)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list, chunksize=CHUNK), total=len(star_descriptions)):
        pass

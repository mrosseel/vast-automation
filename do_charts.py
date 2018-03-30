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
from reading import trash_and_recreate_dir

def set_seaborn_style():
    sns.set_context("notebook", font_scale=4)
    sns.set_style("ticks")

# takes:  [ {'id': star_id, 'match': {'name': match_name, 'separation': separation_deg  } } ]
def plot_lightcurve(tuple, comparison_stars):
#    try:
    star_description = tuple[0]
    star = star_description.local_id
    star_match, separation = get_match(star_description)
    upsilon_text = get_upsilon_string(star_description)
    coord = star_description.coords
    curve = tuple[1]
    if(curve is None):
        print("Curve is None for star", star)
        return
    curve = curve.replace(to_replace=99.99999, value=np.nan, inplace=False)

    curve2_norm = curve
    curve2_norm['V-C'] = curve['V-C'] + comparison_stars[0].vmag

    used_curve = curve2_norm
    used_curve_max = curve2_norm['V-C'].max()
    used_curve_min = curve2_norm['V-C'].min()

    #insert counting column
    used_curve.insert(0, 'Count', range(0, len(used_curve)))
    g = sns.lmplot('Count', 'V-C',
               data=used_curve, size=20, aspect=5,scatter_kws={"s": 15},
               fit_reg=False)
    star_name = '' if star_match == '' else " ({} - dist:{:.4f})".format(star_match, separation)
    plt.title("Star {0}{1}, position: {2}{3}".format(star, star_name, get_hms_dms(coord), upsilon_text))

    #plt.xaxis.set_major_formatter(ticker.FuncFormatter(format_date))
    #plt.set_title("Custom tick formatter")
    #fig.autofmt_xdate()
    plt.xlabel('Count')
    plt.ylabel("Absolute Mag, comp star = {:2.2f}".format(comparison_stars[0].vmag))
    plot_max = used_curve_max
    plot_min = min(plot_max-1, used_curve_min)
    print('min', plot_min, 'max', plot_max, 'usedmin', used_curve_min, 'usedmax', used_curve_max)
    if np.isnan(plot_max) or np.isnan(plot_min):
        print("star is nan", star)
        return
    #print("Star:{},dim:{},bright:{}".format(star, plot_dim, plot_bright))
    plt.ylim(plot_min,plot_max)
    plt.xlim(0, len(used_curve))
    plt.gca().invert_yaxis()
    #g.map(plt.errorbar, 'Count', 'V-C', yerr='s1', fmt='o')
    #plt.ticklabel_format(style='plain', axis='x')
    g.savefig(init.chartsdir+str(star).zfill(5) )
    plt.close(g.fig)
#    except:
#        print("error", tuple)

def plot_phase_diagram(star_description, comparison_stars, suffix='', period=None):
    star = star_description.local_id
    star_match, separation = get_match(star_description)
    upsilon_text = get_upsilon_string(star_description)
    match_string = "({})".format(star_match) if not star_match == '' else ''
    print("Calculating phase diagram for", star)
    curve = reading.read_lightcurve(star)
    if curve is None:
        print("Curve of star {} is None".format(star))
        return
    t_np = curve['JD'].as_matrix()
    y_np = curve['V-C'].as_matrix()
    dy_np = curve['s1'].as_matrix()
    if period is None:
        ls = LombScargleFast()
        period_max = np.max(t_np)-np.min(t_np)
        ls.optimizer.period_range = (0.01,period_max)
        ls.fit(t_np,y_np)
        period = ls.best_period
    print("Best period: " + str(period) + " days")
    fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    plt.xlabel("Phase")
    plt.ylabel("Magnitude")
    plt.title("Star {} {}, period: {:.5f} d{}".format(star, match_string, period, upsilon_text))

    # plotting + calculation of 'double' phase diagram from -1 to 1
    phased_t = np.fmod(t_np/period,1)
    minus_one = lambda t: t - 1
    minus_oner = np.vectorize(minus_one)
    phased_t2 = minus_oner(phased_t)
    phased_lc = y_np[:]
    phased_t_final = np.append(phased_t2, phased_t)
    phased_lc_final = np.append(phased_lc, phased_lc)
    phased_lc_final = phased_lc_final + comparison_stars[0].vmag
    phased_err = np.append(dy_np, dy_np)
    plt.gca().invert_yaxis()
    plt.errorbar(phased_t_final,phased_lc_final,yerr=phased_err,linestyle='none',marker='o', ecolor='gray', elinewidth=1)
    fig.savefig(init.phasedir+'phase'+str(star).zfill(5)+suffix)
    plt.close(fig)

def get_hms_dms(coord):
    return "{:2.0f}$^h$ {:02.0f}$^m$ {:02.2f}$^s$ | {:2.0f}$\degree$ {:02.0f}$'$ {:02.2f}$''$"\
        .format(coord.ra.hms.h, abs(coord.ra.hms.m), abs(coord.ra.hms.s),
                coord.dec.dms.d, abs(coord.dec.dms.m), abs(coord.dec.dms.s))

def format_date(x, pos=None):
    thisind = np.clip(int(x + 0.5), 0, N - 1)
    return r.date[thisind].strftime('%Y-%m-%d')

def store_curve_and_pos(star, star_descriptions, comparison_stars):
    star_description = [x for x in star_descriptions if x.local_id == star][0]
    try:
        tuple = star_description, reading.read_lightcurve(star,filter=False)
        plot_lightcurve(tuple, comparison_stars)
    except FileNotFoundError:
        print("File not found error in store and curve for star", star)

# extract matching strings from star_descr
def get_match(star_description):
    if not star_description.match == None:
        name = star_description.match[0]['catalog_dict']['name']
        separation = star_description.match[0]['separation']
    else:
        name = ''
        separation = ''
    return name, separation

# extract upsilon strings from star_descr
def get_upsilon_string(star_description):
    upsilon = star_description.upsilon
    if not upsilon == None:
        upsilon_text = ", Var: prob={0:.2f},type={1}".format(upsilon['probability'], upsilon['vartype'])
    else:
        upsilon_text = ''
    return upsilon_text


# takes:  [ {'id': star_id, 'match': {'name': match_name, 'separation': separation_deg  } } ]
def run(star_descriptions, comparison_stars):
    star_list = [star.local_id for star in star_descriptions]
    curve_and_pos = []
    set_seaborn_style()
    pool = mp.Pool(init.nr_threads)
    trash_and_recreate_dir(init.chartsdir)

    print("Reading and plotting star positions, total size = ",len(star_list))
    func = partial(store_curve_and_pos, star_descriptions=star_descriptions, comparison_stars=comparison_stars)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list), total=len(star_list)):
        pass

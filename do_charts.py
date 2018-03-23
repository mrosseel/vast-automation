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
def plot_lightcurve(tuple):
#    try:
    star_description = tuple[0]
    star = star_description.local_id
    curve = tuple[1]
    if not star_description.match == None:
        star_match = star_description.match[0]['catalog_dict']['name']
        separation = star_description.match[0]['separation']
    else:
        star_match = ''
        separation = ''

    coord = star_description.coords

    if(curve is None):
        print("Curve is None for star", star)
        return
    curve = curve.replace(to_replace=99.99999, value=np.nan, inplace=False)

    curve_min = curve['V-C'].min()
    curve_max = curve['V-C'].max()
    curve2_norm = curve
    #print("star, min, max:",star, curve_min,curve_max)
    curve2_norm['V-C'] = curve['V-C'] - curve_min

    used_curve = curve2_norm
    used_curve_max = curve2_norm['V-C'].max()

    #insert counting column
    used_curve.insert(0, 'Count', range(0, len(used_curve)))
    g = sns.lmplot('Count', 'V-C',
               data=used_curve, size=20, aspect=5,scatter_kws={"s": 10},
               fit_reg=False)
    plt.title('Star '+ str(star) + ', ' + str(star_match) +', position: ' + get_hms_dms(coord)  +", distance: " + str(separation))

    #plt.xaxis.set_major_formatter(ticker.FuncFormatter(format_date))
    #plt.set_title("Custom tick formatter")
    #fig.autofmt_xdate()
    plt.xlabel('Count')
    plt.ylabel('Mag relative to minimum ' + str(curve_max))
    plt.ylim(max(2, used_curve_max),0)
    plt.xlim(0, len(used_curve))
    plt.gca().invert_yaxis()
    #g.map(plt.errorbar, 'Count', 'V-C', yerr='s1', fmt='o')
    #plt.ticklabel_format(style='plain', axis='x')
    g.savefig(init.chartsdir+str(star).zfill(5) )
    plt.close(g.fig)
#    except:
#        print("error", tuple)

def plot_phase_diagram(star_description, suffix='', period=None):
    star = star_description.local_id
    print("Calculating phase diagram for", star)
    if period is None:
        curve = reading.read_lightcurve(star)
        print(curve)
        if curve is None:
            return
        t_np = curve['JD'].as_matrix()
        y_np = curve['V-C'].as_matrix()
        dy_np = curve['s1'].as_matrix()
        ls = LombScargleFast()
        period_max = np.max(t_np)-np.min(t_np)
        ls.optimizer.period_range = (0.01,period_max)
        ls.fit(t_np,y_np)
        period = ls.best_period
    print("Best period: " + str(period) + " days")
    fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    plt.xlabel("Phase")
    plt.ylabel("Diff mag")
    plt.title("Lomb-Scargle-Periodogram for star "+str(star)+" period: " + str(period))

    # plotting + calculation of 'double' phase diagram from -1 to 1
    phased_t = np.fmod(t_np/period,1)
    minus_one = lambda t: t - 1
    minus_oner = np.vectorize(minus_one)
    phased_t2 = minus_oner(phased_t)
    phased_lc = y_np[:]
    phased_t_final = np.append(phased_t2, phased_t)
    phased_lc_final = np.append(phased_lc, phased_lc)
    phased_err = np.append(dy_np, dy_np)
    plt.errorbar(phased_t_final,phased_lc_final,yerr=phased_err,linestyle='none',marker='o')
    fig.savefig(init.phasedir+'phase'+str(star).zfill(5)+suffix)
    plt.close(fig)

def get_hms_dms(coord):
    return str(coord.ra.hms.h) + ' ' + str(abs(coord.ra.hms.m)) + ' ' + str(abs(round(coord.ra.hms.s, 2))) \
           + ' | ' + str(coord.dec.dms.d) + ' ' + str(abs(coord.dec.dms.m)) + ' ' + str(abs(round(coord.dec.dms.s, 1)))

def format_date(x, pos=None):
    thisind = np.clip(int(x + 0.5), 0, N - 1)
    return r.date[thisind].strftime('%Y-%m-%d')

def store_curve_and_pos(star, star_descriptions):
    star_description = [x for x in star_descriptions if x.local_id == star][0]
    try:
        tuple = star_description, reading.read_lightcurve(star,filter=False)
        plot_lightcurve(tuple)
    except FileNotFoundError:
        print("File not found error in store and curve for star", star)

# takes:  [ {'id': star_id, 'match': {'name': match_name, 'separation': separation_deg  } } ]
def run(star_descriptions):
    star_list = [star.local_id for star in star_descriptions]
    curve_and_pos = []
    set_seaborn_style()
    pool = mp.Pool(init.nr_threads)
    trash_and_recreate_dir(init.chartsdir)

    print("Reading and plotting star positions, total size = ",len(star_list))
    func = partial(store_curve_and_pos, star_descriptions=star_descriptions)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list), total=len(star_list)):
        pass

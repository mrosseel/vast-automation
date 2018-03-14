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
import init
from reading import trash_and_recreate_dir

def set_seaborn_style():
    sns.set_context("notebook", font_scale=1.1)
    sns.set_style("ticks")

def plot_lightcurve(tuple, matches):
#    try:
    star = tuple[0]
    curve = tuple[1]
    pos = tuple[2]
    star_match = matches[star][4]
    separation = matches[star][7]
    coord = SkyCoord(pos[0], pos[1], unit='deg')

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

    #insert counting column
    used_curve.insert(0, 'Count', range(0, len(used_curve)))
    g = sns.lmplot('Count', 'V-C',
               data=used_curve, size=20, aspect=5,scatter_kws={"s": 10},
               fit_reg=False)
    #print(used_curve.head(10))
    #print(coord.ra.hms, coord.dec.dms)
    plt.title('Star '+ str(star) + " : " + str(coord.ra.hms) + ' - ' + str(coord.dec.dms) + " - " + str(star_match) + ", " + str(separation))

    #plt.xaxis.set_major_formatter(ticker.FuncFormatter(format_date))
    #plt.set_title("Custom tick formatter")
    #fig.autofmt_xdate()
    plt.xlabel('Count')
    plt.ylabel('Mag relative to minimum ' + str(curve_max))
    plt.ylim(2,0)
    plt.xlim(0, len(used_curve))
    plt.gca().invert_yaxis()
    #g.map(plt.errorbar, used_curve['Count'], used_curve['V-C'], yerr=used_curve['s1'], fmt='o')
    #plt.ticklabel_format(style='plain', axis='x')
    g.savefig(init.chartsdir+str(star).zfill(5) )
    plt.close(g.fig)
#    except:
#        print("error", tuple)

def format_date(x, pos=None):
    thisind = np.clip(int(x + 0.5), 0, N - 1)
    return r.date[thisind].strftime('%Y-%m-%d')

def store_curve_and_pos(star):
    try:
        tuple = star, reading.read_lightcurve(star,filter=False), reading.read_worldpos(star)
        return tuple
    except FileNotFoundError:
        print("File not found error in store and curve for star", star)

def run(matches):
    star_list = [*matches] # unpack
    curve_and_pos = []
    set_seaborn_style()
    pool = mp.Pool(init.nr_threads)
    print("Reading star positions, total size = ",len(star_list))
    for _ in tqdm.tqdm(pool.imap_unordered(store_curve_and_pos, star_list), total=len(init.all_star_list)):
        curve_and_pos.append(_)
        pass
    print("Plotting stars, total size = ",len(curve_and_pos))
    trash_and_recreate_dir(init.chartsdir)
    func = partial(plot_lightcurve, matches=matches)
    for _ in tqdm.tqdm(pool.imap_unordered(func, curve_and_pos), total=len(curve_and_pos)):
        pass

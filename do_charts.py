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
    sns.set_context("notebook", font_scale=4)
    sns.set_style("ticks")

# takes:  [ {'id': star_id, 'match': {'name': match_name, 'separation': separation_deg  } } ]
def plot_lightcurve(tuple):
#    try:
    star = tuple[0]
    curve = tuple[1]
    pos = tuple[2]
    if len(tuple) == 4:
        star_match = tuple[3]['name']
        separation = tuple[3]['separation']
    else:
        star_match = ''
        separation = ''

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

def get_hms_dms(coord):
    return str(coord.ra.hms.h) + ' ' + str(abs(coord.ra.hms.m)) + ' ' + str(abs(round(coord.ra.hms.s, 2))) \
           + ' | ' + str(coord.dec.dms.d) + ' ' + str(abs(coord.dec.dms.m)) + ' ' + str(abs(round(coord.dec.dms.s, 1)))

def format_date(x, pos=None):
    thisind = np.clip(int(x + 0.5), 0, N - 1)
    return r.date[thisind].strftime('%Y-%m-%d')

def store_curve_and_pos(chart_object):
    try:
        print(chart_object)
        star = chart_object['id']
        tuple = star, reading.read_lightcurve(star,filter=False), reading.read_worldpos(star)
        if 'match' in chart_object.keys():
            tuple= tuple + (chart_object['match'],)
        return tuple
    except FileNotFoundError:
        print("File not found error in store and curve for star", star)

# takes:  [ {'id': star_id, 'match': {'name': match_name, 'separation': separation_deg  } } ]
def run(matches):
    star_list = [*matches] # unpack
    curve_and_pos = []
    set_seaborn_style()
    pool = mp.Pool(init.nr_threads)
    print("Reading star positions, total size = ",len(star_list))

    func = partial(store_curve_and_pos)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list), total=len(star_list)):
        curve_and_pos.append(_)
        pass

    print("Plotting stars, total size = ",len(curve_and_pos))
    trash_and_recreate_dir(init.chartsdir)
    func = partial(plot_lightcurve)
    for _ in tqdm.tqdm(pool.imap_unordered(func, curve_and_pos), total=len(curve_and_pos)):
        pass

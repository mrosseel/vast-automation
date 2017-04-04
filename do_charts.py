import init
import reading
import matplotlib as mp
mp.use('Agg') # needs no X server
import matplotlib.pyplot as plt
import seaborn as sns
import os
import multiprocessing as mp
import tqdm


def set_seaborn_style():
    sns.set_context("notebook", font_scale=1.1)
    sns.set_style("ticks")

def plot_lightcurve(tuple):
    try:
        star = tuple[0]
        curve = tuple[1]
        pos = tuple[2]

        curve_min = curve['V-C'].min()
        curve_max = curve['V-C'].max()
        curve2 = curve
        #print("min, max:",curve_min,curve_max)
        curve2['V-C'] = curve['V-C'] - curve_min

        #insert counting column
        curve2.insert(0, 'Count', range(0, len(curve2)))
        g = sns.lmplot('Count', 'V-C',
                   data=curve2,
                   fit_reg=False)

        plt.title('Star '+ str(star))
        #+ " : " + pixel_to_radec(wcs_config, pos[0], pos[1]).to_string('hmsdms') + ' - ' +str(pos[0]) + ', ' + str(pos[1]))
        plt.xlabel('Obs #')
        plt.ylabel('Mag')
        plt.ylim(2,0)
        plt.gca().invert_yaxis()
        #plt.ticklabel_format(style='plain', axis='x')
        #fig = plt.figure(figsize=(70,10))
        #sns.plt.show()
        g.savefig(init.lightcurve_dir+str(star).zfill(5) )
        plt.close(g.fig)
    except:
        print("error", tuple)


def store_curve_and_pos(star):
    try:
        tuple = star, reading.read_lightcurve(star), reading.read_pos(star)
        return tuple
    except FileNotFoundError:
        print("File not found error in store and curve for star", star)


def run():
    curve_and_pos = []
    set_seaborn_style()
    pool = mp.Pool(8)
    star_list = init.all_star_list
    print("Reading star positions, total size = ",len(star_list))
    for _ in tqdm.tqdm(pool.imap_unordered(store_curve_and_pos, star_list), total=len(init.all_star_list)):
        curve_and_pos.append(_)
        pass
    print("Plotting stars, total size = ",len(curve_and_pos))
    for _ in tqdm.tqdm(pool.imap_unordered(plot_lightcurve, curve_and_pos), total=len(curve_and_pos)):
        pass

run()

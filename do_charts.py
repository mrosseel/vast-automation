import init
import pandas as pd
import numpy as np
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

def read_lightcurve(star):
    #print("Reading lightcurve", star)
    df = pd.read_csv(init.lightcurve_dir + 'curve_' + str(star).zfill(5) + '.txt', skiprows=[1], sep=' ')
    df = df[df['V-C'] < 99]
    return preprocess_lightcurve(df)

def preprocess_lightcurve(df):
    try:
        P = np.percentile(df['V-C'], [5, 95])
        df2 = df[(df['V-C'] > P[0]) & (df['V-C'] < P[1])]
        return df2
    except IndexError:
        print("len df:", len(df))

def read_pos(star):
    try:
        df = pd.read_csv(init.lightcurve_dir + 'pos_' + str(star).zfill(5) + '.txt', skiprows=[1], sep=' ')
        df2 = df[df['X'] > 0]
        df3 = df2[df['MAG'] < 99]
        return (df3['X'].iloc[0], df3['Y'].iloc[0])
    except IndexError:
        print("df:",len(df),"df2:", len(df2),"df3:", len(df3))
        print(len(df))

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
        plt.ticklabel_format(style='plain', axis='x')
        #fig = plt.figure(figsize=(70,10))
        #sns.plt.show()
        g.savefig(init.lightcurve_dir+str(star).zfill(5) )
        plt.close(g.fig)
    except:
        print("error", tuple)


def store_curve_and_pos(star):
    try:
        tuple = star, read_lightcurve(star), read_pos(star)
        return tuple
    except FileNotFoundError:
        print("File not found error in store and curve for star", star)


def run():
    curve_and_pos = []

    set_seaborn_style()
    pool = mp.Pool(8)
    #init.all_star_list = range(1,10)
    print("Reading star positions, total size = ",len(init.all_star_list))
    for _ in tqdm.tqdm(pool.imap_unordered(store_curve_and_pos, init.all_star_list), total=len(init.all_star_list)):
        curve_and_pos.append(_)
        pass
    print("Plotting stars, total size = ",len(curve_and_pos))
    for _ in tqdm.tqdm(pool.imap_unordered(plot_lightcurve, curve_and_pos), total=len(curve_and_pos)):
        pass

run()

#try:
#    for star in init.all_star_list:
#        plot_lightcurve(read_lightcurve(star), read_pos(star), star)
#except D:
#    print("D", D)

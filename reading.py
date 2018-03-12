import init
import os
import shutil
import pandas as pd
import numpy as np
import sys
import errno

def read_lightcurve(star,filter=True,preprocess=True):
    try:
        #print("Reading lightcurve", star, init.lightcurve_dir + 'curve_' + str(star).zfill(5) + '.txt')
        df = pd.read_csv(init.lightcurvedir + 'curve_' + str(star).zfill(5) + '.txt', skiprows=[1], sep=' ')
        if(filter):
            df = df[df['V-C'] < 99]
        if(preprocess):
            df = preprocess_lightcurve(df)
        return df
    except OSError:
        print("OSError for star:", star)

def preprocess_lightcurve(df):
    try:
        P = np.percentile(df['V-C'], [5, 95])
        df2 = df[(df['V-C'] > P[0]) & (df['V-C'] < P[1])]
        return df2
    except IndexError:
        print("len df:", len(df))



def read_pos_old(star):
    return "TODO: position is not yet returned"
    try:
        df = pd.read_csv(init.posdir + 'pos_' + str(star).zfill(5) + '.txt', skiprows=[1], sep=' ')
        df2 = df[df['X'] > 0]
        df3 = df2[df['MAG'] < 99]
        return (df3['X'].iloc[0], df3['Y'].iloc[0])
    except IndexError:
        print("ERROR: IndexError")
        #print("df:",len(df),"df2:", len(df2),"df3:", len(df3))
        print(len(df))

def read_pos(star, jd):
    try:
        df = pd.read_csv(init.posdir + 'pos_' + str(star).zfill(5) + '.txt', skiprows=[1], sep=' ')
        print(df.head())
        df2 = df[df['X'] > 0]
        df3 = df2[df['MAG'] < 99]
        row = df.loc[df['JD'] == jd]
        print("row", row, jd)
        row = df3.iloc[0]
        return [row['JD'], row['X'],row['Y'], row['MAG']]
        #return (df3['X'].iloc[0], df3['Y'].iloc[0])
        return df
    except IndexError:
        print("ERROR: IndexError")
        #print("df:",len(df),"df2:", len(df2),"df3:", len(df3))
        print(len(df))

def trash_and_recreate_dir(dir):
    shutil.rmtree(dir, ignore_errors=True)
    os.makedirs(dir, exist_ok=True)

# searches for the last written star in the path, and returns a star list including that star so it can be overwritten
def search_last_star(star_list, the_path):
    the_dir = os.listdir(the_path)
    the_dir.sort()
    import re
    m = re.search('\d+',the_dir[-1])
    last_star = int(m.group(0).lstrip('0'))
    last_star_index = star_list.index(last_star)
    return star_list[last_star_index:]

# read the world positions and return them in a nicely formatted list
def read_world_positions(the_path):
    the_dir = os.listdir(the_path)
    the_dir.sort()
    results = {}
    for name in the_dir: # 'file' is a builtin type, 'name' is a less-ambiguous variable name.
        try:
            with open(the_path + name) as f: # No need to specify 'r': this is the default.
                results[name] = f.readlines()[0].split(' ')
        except IOError as exc:
            if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
                raise # Propagate other kinds of IOError.
    return results


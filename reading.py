import init
import os
import shutil
import pandas as pd
import numpy as np
import errno
import re


def read_lightcurve(star,filter=True,preprocess=True, directory=init.lightcurvedir):
    try:
        #print("Reading lightcurve", star, init.lightcurve_dir + 'curve_' + str(star).zfill(5) + '.txt')
        df = pd.read_csv(directory + 'curve_' + str(star).zfill(5) + '.txt', skiprows=[1], sep=' ')
        if(filter):
            df = df[df['V-C'] < 99]
        if(preprocess):
            df = preprocess_lightcurve(df)
        return df
    except OSError as e:
        print("OSError for star:", star, e)

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

def read_worldpos(star):
    with open(get_worldpos_filename(star)) as f: # No need to specify 'r': this is the default.
        return f.readlines()[0].split(' ')

def trash_and_recreate_dir(dir):
    os.system('rm -fr "%s"' % dir)
    #shutil.rmtree(dir, ignore_errors=True)
    os.makedirs(dir, exist_ok=True)

# helper function
def get_pos_filename(star):
    return init.posdir + "pos_" + str(star).zfill(5) + ".txt"

def get_worldpos_filename(star):
    return init.worldposdir + "worldpos_" + str(star).zfill(5) + ".txt"

# takes a star_list and a dir, and returns a reduced star list - all stars which already have a file in that dir are removed
def reduce_star_list(star_list, the_path):
    the_dir = os.listdir(the_path)
    the_dir.sort()
    found = []
    for filename in the_dir:
        found.append(filename_to_star(filename))
    print("Found", len(found), "stars already processed in", the_path)
    return [item for item in star_list if item not in found]

# takes a filename and extracts the star number from it
def filename_to_star(filename):
    import re
    m = re.search('\d+',filename)
    return int(m.group(0).lstrip('0'))

# read the world positions and return them in a dictionary
# returns {'name': [ra.deg, dec.deg ]}
def read_world_positions(the_path):
    the_dir = os.listdir(the_path)
    the_dir.sort()
    results = {}
    for name in the_dir: # 'file' is a builtin type, 'name' is a less-ambiguous variable name.
        try:
            with open(the_path + name) as f: # No need to specify 'r': this is the default.
                results[filename_to_star(name)] = f.readlines()[0].split(' ')
        except IOError as exc:
            if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
                raise # Propagate other kinds of IOError.
    return results

def read_comparison_star():
    with open(init.basedir + 'munifind.txt', 'r') as fp:
        for i, line in enumerate(fp):
            if i == 1:
                print(line)# 26th line
                m = re.search('Reference star:\s*(\d+),', line)
                comparison_star = int(m.group(1))
                break
    return comparison_star

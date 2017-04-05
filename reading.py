import init
import pandas as pd
import numpy as np

def read_lightcurve(star,filter=True,preprocess=True):
    #print("Reading lightcurve", star, init.lightcurve_dir + 'curve_' + str(star).zfill(5) + '.txt')
    df = pd.read_csv(init.lightcurve_dir + 'curve_' + str(star).zfill(5) + '.txt', skiprows=[1], sep=' ')
    if(filter):
        df = df[df['V-C'] < 99]
    if(preprocess):
        df = preprocess_lightcurve(df)
    return df

def preprocess_lightcurve(df):
    try:
        P = np.percentile(df['V-C'], [5, 95])
        df2 = df[(df['V-C'] > P[0]) & (df['V-C'] < P[1])]
        return df2
    except IndexError:
        print("len df:", len(df))



def read_pos(star):
    return "TODO: position is not yet returned"
    try:
        df = pd.read_csv(init.lightcurve_dir + 'pos_' + str(star).zfill(5) + '.txt', skiprows=[1], sep=' ')
        df2 = df[df['X'] > 0]
        df3 = df2[df['MAG'] < 99]
        return (df3['X'].iloc[0], df3['Y'].iloc[0])
    except IndexError:
        print("ERROR: IndexError")
        #print("df:",len(df),"df2:", len(df2),"df3:", len(df3))
        print(len(df))

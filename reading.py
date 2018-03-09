import init
import pandas as pd
import numpy as np

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

read_pos(1, 1)

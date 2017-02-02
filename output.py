import init
import astropy_helper
import upsilon_helper
import pandas as pd
import os
import numpy as np

# import matplotlib as mp
# import matplotlib.pyplot as plt
# import seaborn as sns

#cell 14
import subprocess
from subprocess import call
# munilist -a 2 --object 143 -v 143 -

select_star_list = [143,6394,598,2675,3111,2584]
single_star_list = [598]
all_star_list = range(1,1000)
star_list = all_star_list




#cell 3
#!pwd
#!rm {init.basedir+'*.fts'}
#!konve {init.basedir+'*.fit'} -o {init.basedir+'kout??????.fts'}
#!muniphot {init.basedir+'*.fts'} -p muniphot.conf -o {init.basedir+'phot??????.pht'}
#!munimatch -s sp_fields=1 {init.basedir+'phot0000001.pht'} {init.basedir+'phot??????.pht'} -o {init.basedir+'match??????.pht'}
#!munifind -a {aperture} {init.basedir+'munifind.txt'} {init.basedir+'match*'}

# INDEX MEAN_MAG STDEV GOODPOINTS
def read_munifind(filename):
    df = pd.read_csv(filename, skiprows=[1], sep=' ')
    df.rename(columns = {'INDEX':'STAR'}, inplace = True)
    print("max goodpoints:", df['GOODPOINTS'].max())
    print("min stdev:", df['STDEV'].min())
    print(df.sort_values('STDEV').head())
    return df

#cell 11
def getBestComparisonStars(df):
    result = []
    # TODO: filter all stars which have not the maximum of GOODPOINTS
    df_lowest_stdev = df.sort_values('STDEV')
    return df_lowest_stdev.head(10)

def do_best_comps(df):
    bestcomps = getBestComparisonStars(df)
    print(bestcomps)
    check_stars = []
    for index, row in bestcomps.iterrows():
        #print(row, '\n')
        check_stars.append(int(row['STAR']))
    check_stars_str = ','.join(map(str, check_stars))
    print(check_stars_str)

    #cell 12
    print(df.iloc[143])

    return check_stars_str

def write_lightcurve(checkstar, star_list, aperture, lightcurve_dir, check_stars_str):
    for star in star_list:
        print("lightcurve:", star)
#        print("--verbose -a ", str(aperture), " -q --object ", str(star), " -v ", str(star),
#              " -c ", str(check_stars_str), (lightcurve_dir + str(star) + ".txt"), (init.basedir+'match*.pht'))
        os.system("munilist -a "+str(aperture)+ " -q --object "+ str(star)+ " -v "+ str(star)+ " -c "+ check_stars_str+ " " + init.lightcurve_dir + str(star) + ".txt "+ init.basedir+'match*.pht >/dev/null')
#        !munilist --verbose -a {str(aperture)} -q --object {str(star)} -v {str(star)} -c {str(8)} {lightcurve_dir + str(star) + ".txt"} {init.basedir+'match*.pht'}

def write_pos(star_list, aperture, lightcurve_dir):
    for star in star_list:
        print("pos:", star)
        os.system("munilist -a " + str(aperture)+ " -q --obj-plot --object "+ str(star)+ " " + init.lightcurve_dir + "pos_" + str(star) + ".txt "+ init.basedir+'match*.pht >/dev/null')

def do_write_post_and_curve(df, check_stars_str):
    #star_list = (143,264,2675,1045,847,1193)
    call(["mkdir", init.lightcurve_dir])
    write_pos(star_list, init.aperture, init.lightcurve_dir)
    write_lightcurve(6, star_list, init.aperture, init.lightcurve_dir, check_stars_str) # TODO hard-coded checkstar


df = read_munifind(init.basedir+'munifind.txt')
check_stars_str = do_best_comps(df)
do_write_post_and_curve(df, check_stars_str)
upsilon_helper.predict_star_list(all_star_list)

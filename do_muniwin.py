import init
import astropy_helper
import upsilon_helper
import pandas as pd
import os
import numpy as np
import multiprocessing as mp
import tqdm

# import matplotlib as mp
# import matplotlib.pyplot as plt
# import seaborn as sns

#cell 14
import subprocess
from subprocess import call
# munilist -a 2 --object 143 -v 143 -

#select_star_list = [143,6394,598,2675,3111,2584]
#single_star_list = [598]
#all_star_list = range(1,10000)
#star_list = all_star_list

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
    maxPoints = df['GOODPOINTS'].max()
    df = df[df['GOODPOINTS'] > maxPoints*0.99]
    df_lowest_stdev = df.sort_values('STDEV')
    comparison_stars = df_lowest_stdev.head(10)
    print("Comparison stars: ", comparison_stars)
    return comparison_stars

def do_best_comps(df):
    bestcomps = getBestComparisonStars(df)
    print(bestcomps)
    check_stars = []
    for index, row in bestcomps.iterrows():
        #print(row, '\n')
        check_stars.append(int(row['STAR']))
    check_stars_str = ','.join(map(str, check_stars))
    print(check_stars_str)
    return check_stars_str

def write_photometry():
    print("write photometry")
    # !rm {init.basedir+'*.fts'}
    os.system('rm '+init.basedir+'*.fts')
    # !konve {init.basedir+'*.fit'} -o {init.basedir+'kout??????.fts'}
    os.system('konve '+init.basedir+'*.fit -o '+init.basedir+'kout??????.fts')
    # !muniphot {init.basedir+'*.fts'} -p muniphot.conf -o {init.basedir+'phot??????.pht'}
    os.system('muniphot '+init.basedir+'*.fts -p muniphot.conf -o '+init.basedir+'phot??????.pht')

def write_match(base_photometry_file):
    print("write match", base_photometry_file)
    os.system('rm '+init.basedir+'match*')
    # !munimatch -s sp_fields=1 {init.basedir+'phot0000001.pht'} {init.basedir+'phot??????.pht'} -o {init.basedir+'match??????.pht'}
    os.system('munimatch -s sp_fields=1 '+base_photometry_file+' '+init.basedir+'phot??????.pht -o '+init.basedir+'match??????.pht')

def write_munifind():
    print("write munifind")
    # !munifind -a {aperture} {init.basedir+'munifind.txt'} {init.basedir+'match*'}
    os.system('munifind -a '+str(init.aperture)+' '+init.basedir+'munifind.txt '+init.basedir+'match*')

def write_munifind_check_stars():
    print("write munifind")
    # !munifind -a {aperture} {init.basedir+'munifind.txt'} {init.basedir+'match*'}
    os.system('munifind -a '+str(init.aperture)+' '+' -c '+ check_stars_str+' '+init.basedir+'munifind.txt '+init.basedir+'match*')



def write_lightcurve(star):
#        print("--verbose -a ", str(aperture), " -q --object ", str(star), " -v ", str(star),
#              " -c ", str(check_stars_str), (lightcurve_dir + str(star) + ".txt"), (init.basedir+'match*.pht'))
    os.system("munilist -a "+str(init.aperture)+ " -q --object "+ str(star)+ " -v "+ str(star)+ " -c "+ check_stars_str+ " " + init.lightcurve_dir + str(star) + ".txt "+ init.basedir+'match*.pht >/dev/null')
#        !munilist --verbose -a {str(aperture)} -q --object {str(star)} -v {str(star)} -c {str(8)} {lightcurve_dir + str(star) + ".txt"} {init.basedir+'match*.pht'}

def write_pos(star):
    os.system("munilist -a " + str(init.aperture)+ " -q --obj-plot --object "+ str(star)+ " " + init.lightcurve_dir + "pos_" + str(star) + ".txt "+ init.basedir+'match*.pht >/dev/null')

def do_write_post_and_curve(df, star_list, check_stars_str):
    call(["mkdir", init.lightcurve_dir])
    pool = mp.Pool(8)
    print("Writing star positions for ",len(star_list),"stars into ",init.lightcurve_dir)
    for _ in tqdm.tqdm(pool.imap_unordered(write_pos, star_list), total=len(star_list)):
        pass
    for _ in tqdm.tqdm(pool.imap_unordered(write_lightcurve, star_list), total=len(star_list)):
        pass


#write_photometry()
#write_match(init.basedir+'phot000046.pht')
write_munifind()
df = read_munifind(init.basedir+'munifind.txt')
check_stars_str = do_best_comps(df)
write_munifind_check_stars()
df = read_munifind(init.basedir+'munifind.txt')
#star_list = (143,264,2675,1045,847,1193)
star_list = init.all_star_list
do_write_post_and_curve(df, init.all_star_list, check_stars_str)

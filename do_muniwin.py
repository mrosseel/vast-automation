import init
import astropy_helper
#import upsilon_helper
import pandas as pd
import numpy as np
import multiprocessing as mp
import logging
import tqdm
import os, sys
import subprocess
from subprocess import call
from functools import partial

# INDEX MEAN_MAG STDEV GOODPOINTS
def read_munifind(filename):
    df = pd.read_csv(filename, skiprows=[1], sep=' ')
    df.rename(columns = {'INDEX':'STAR'}, inplace = True)
    print("max goodpoints:", df['GOODPOINTS'].max())
    print("min stdev:", df['STDEV'].min())
    print(df.sort_values('STDEV').head())
    return df

def getBestComparisonStars(df):
    result = []
    maxPoints = df['GOODPOINTS'].max()
    df = df[df['GOODPOINTS'] > maxPoints*0.99]
    df_lowest_stdev = df.sort_values('STDEV')
    comparison_stars = df_lowest_stdev.head(12)
    print("Comparison stars: ", comparison_stars)
    return comparison_stars

def do_best_comps(bestcomps):
    check_stars = []
    for index, row in bestcomps.iterrows():
        #print(row, '\n')
        check_stars.append(int(row['STAR']))
    return check_stars

def join_check_stars(check_stars, exclude_star):
    check_stars = filter(lambda star: star != exclude_star, check_stars)
    check_stars_string = ','.join(map(str, check_stars))
    return check_stars_string

def write_convert_fits():
    print("convert fits files to fts")
    # !rm {init.basedir+'*.fts'}
    os.system('rm '+init.basedir+'*.fts')
    # !konve {init.basedir+'*.fit'} -o {init.basedir+'kout??????.fts'}
    os.system('konve '+init.basedir+'*.fit -o '+init.basedir+'kout??????.fts')

def write_photometry():
    print("write photometry")
    os.system('rm '+init.basedir+'*.pht')
    # !muniphot {init.basedir+'*.fts'} -p muniphot.conf -o {init.basedir+'phot??????.pht'}
    os.system('muniphot '+init.basedir+'*.fts -p muniphot.conf -o '+init.basedir+'phot??????.pht')

def write_match(base_photometry_file):
    print("write match with base file:", base_photometry_file)
    os.system('rm '+init.basedir+'match*')
    # !munimatch -s sp_fields=1 {init.basedir+'phot0000001.pht'} {init.basedir+'phot??????.pht'} -o {init.basedir+'match??????.pht'}
    os.system('munimatch -s sp_fields=1 '+base_photometry_file+' '+init.basedir+'phot??????.pht -o '+init.basedir+'match??????.pht')

def write_munifind():
    print("write munifind")
    # !munifind -a {aperture} {init.basedir+'munifind.txt'} {init.basedir+'match*'}
    os.system('munifind -a '+str(init.aperture)+' '+init.basedir+'munifind.txt '+init.basedir+'match*')

def write_munifind_check_stars(check_star):
    print("write munifind check stars using check star:", check_star)
    # !munifind -a {aperture} {init.basedir+'munifind.txt'} {init.basedir+'match*'}
    os.system('munifind -a '+str(init.aperture)+' '+' -c '+ str(check_star) +' '+init.basedir+'munifind.txt '+init.basedir+'match*')



def write_lightcurve(star, check_stars_list):
    check_stars = join_check_stars(check_stars_list, star)
    os.system("munilist -a "+str(init.aperture)+ " -q --object "+ str(star)+ " -v "+ str(star)+ " -c "+ check_stars + " " + init.lightcurve_dir + 'curve_' + str(star).zfill(5) + ".txt "+ init.basedir+'match*.pht >/dev/null')
    #print("--verbose -a ", str(init.aperture), " -q --object ", str(star), " -v ", str(star), " -c ", check_stars, (init.lightcurve_dir + 'curve_' + str(star).zfill(5) + ".txt "), (init.basedir+'match*.pht  >/dev/null'))
    # !munilist --verbose -a {str(aperture)} -q --object {str(star)} -v {str(star)} -c {str(8)} {lightcurve_dir + str(star) + ".txt"} {init.basedir+'match*.pht'}


#TODO add check stars to this command?
def write_pos(star, check_stars_list):
    check_stars = join_check_stars(check_stars_list, star)
    os.system("munilist -a " + str(init.aperture)+ " -q --obj-plot --object "+ str(star)+ " " + init.lightcurve_dir + "pos_" + str(star).zfill(5) + ".txt "+ init.basedir+'match*.pht >/dev/null')

def do_write_pos(star_list, check_stars_list):
    call(["mkdir", init.lightcurve_dir])
    pool = mp.Pool(8)
    func = partial(write_pos, check_stars_list=check_stars_list)
    print("Writing star positions for ",len(star_list),"stars into ",init.lightcurve_dir)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list), total=len(star_list)):
        pass

def do_write_curve(star_list, check_stars_list):
    call(["mkdir", init.lightcurve_dir])
    pool = mp.Pool(8)
    func = partial(write_lightcurve, check_stars_list=check_stars_list)
    print("Writing star lightcurves for ",len(star_list),"stars into ",init.lightcurve_dir)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list), total=len(star_list)):
        pass

def run_determine_reference_frame():
    write_convert_fits()
    write_photometry()

def run_do_rest(reference_phot):
    write_match(reference_phot)
    write_munifind()
    df = read_munifind(init.basedir+'munifind.txt')
    bestcomps = getBestComparisonStars(df)
    print("bestcomps: ", bestcomps)
    check_stars_list = do_best_comps(bestcomps)
    print("check_stars_list: ", check_stars_list)
    write_munifind_check_stars(check_stars_list[0])
    #star_list = (143,264,2675,1045,847,1193)
    star_list = init.all_star_list
    do_write_curve(star_list, check_stars_list)
    do_write_pos(star_list, check_stars_list)

#logger = mp.log_to_stderr()
#logger.setLevel(mp.SUBDEBUG)
#run_determine_reference_frame()
run_do_rest(init.basedir+'phot000001.pht')

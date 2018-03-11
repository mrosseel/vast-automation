import init
import do_calibration
import do_charts
from init import trash_and_recreate_dir
import pandas as pd
import multiprocessing as mp
import logging
import tqdm
import os, sys
from functools import partial
from subprocess import call
import pickle
import time

# INDEX MEAN_MAG STDEV GOODPOINTS
def read_munifind(filename):
    df = pd.read_csv(filename, skiprows=[1], sep=' ')
    df.rename(columns = {'INDEX':'STAR'}, inplace = True)
    print("max goodpoints:", df['GOODPOINTS'].max())
    print("min stdev:", df['STDEV'].min())
    print(df.sort_values('STDEV').head())
    return df

# gets the stars with a maximum of measurements and lowest stdev
def getBestComparisonStars(nrOfStars):
    df = read_munifind(init.basedir+'munifind_temp.txt')
    maxPoints = df['GOODPOINTS'].max()
    df = df[df['GOODPOINTS'] > maxPoints*0.99]
    df_lowest_stdev = df.sort_values('STDEV')
    comparison_stars = df_lowest_stdev.head(nrOfStars)
    print("Comparison stars: ", comparison_stars)
    return comparison_stars

def do_best_comparison_stars(nrOfStars):
    bestcomps = getBestComparisonStars(nrOfStars)
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
    trash_and_recreate_dir(init.convfitsdir)
    os.system('konve '+init.fitsdir+'*.fit -o '+init.convfitsdir+'kout??????.fts')

def write_photometry():
    print("write photometry")
    trash_and_recreate_dir(init.photometrydir)
    os.system('muniphot '+init.convfitsdir+'*.fts -p muniphot.conf -o '+init.photometrydir+'phot??????.pht')

def write_match(base_photometry_file):
    print("write match with base file:", base_photometry_file)
    trash_and_recreate_dir(init.matchedphotometrydir)
    os.system('munimatch -s sp_fields=1 '+base_photometry_file+' '+init.photometrydir+'phot??????.pht -o '+init.matchedphotometrydir+'match??????.pht')

def write_munifind():
    print("write munifind")
    os.system('munifind -a '+str(init.aperture)+' '+init.basedir+'munifind_temp.txt '+init.matchedphotometrydir+'match*')

def write_munifind_check_stars(check_star):
    print("write munifind check stars using check star:", check_star)
    os.system('munifind -a '+str(init.aperture)+' '+' -c '+ str(check_star) +' '+init.basedir+'munifind.txt '+init.matchedphotometrydir+'match*')

def write_lightcurve(star, check_stars_list):
    check_stars = join_check_stars(check_stars_list, star)
    os.system("munilist -a "+str(init.aperture)+ " -q --object "+ str(star)+ " -v "+ str(star)+ " -c "+ check_stars + " " + init.lightcurvedir + 'curve_' + str(star).zfill(5) + ".txt "+ init.matchedphotometrydir+'match*.pht >/dev/null')
    #print("--verbose -a ", str(init.aperture), " -q --object ", str(star), " -v ", str(star), " -c ", check_stars, (init.lightcurve_dir + 'curve_' + str(star).zfill(5) + ".txt "), (init.basedir+'match*.pht  >/dev/null'))
    # !munilist --verbose -a {str(aperture)} -q --object {str(star)} -v {str(star)} -c {str(8)} {lightcurve_dir + str(star) + ".txt"} {init.basedir+'match*.pht'}


#TODO add check stars to this command?
def write_pos(star, check_stars_list):
    start = time.time()
    check_stars = join_check_stars(check_stars_list, star)
    #os.system("munilist -a " + str(init.aperture)+ " -q --obj-plot --object "+ str(star)+ " " + get_pos_filename(star) + " " + init.matchedphotometrydir+'match*.pht >/dev/null')
    call("munilist -a " + str(init.aperture)+ " -q --obj-plot --object "+ str(star)+ " " + get_pos_filename(star) + " " + init.matchedphotometrydir+'match*.pht >/dev/null', shell=True)
    end = time.time()

def do_write_pos(star_list, check_stars_list):
    trash_and_recreate_dir(init.posdir)
    pool = mp.Pool(8, maxtasksperchild=100)
    func = partial(write_pos, check_stars_list=check_stars_list)
    print("Writing star positions for",len(star_list),"stars into",init.posdir)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list, 10), total=len(star_list)):
        pass

def do_write_curve(star_list, check_stars_list):
    trash_and_recreate_dir(init.lightcurvedir)
    pool = mp.Pool(8, maxtasksperchild=100)
    func = partial(write_lightcurve, check_stars_list=check_stars_list)
    print("Writing star lightcurves for",len(star_list),"stars into",init.lightcurvedir)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list, 10), total=len(star_list)):
        pass

def do_world_pos(wcs, star_list, reference_frame_index):
    trash_and_recreate_dir(init.worldposdir)
    print("index", reference_frame_index)
    pool = mp.Pool(1)
    func = partial(world_pos, wcs=wcs, reference_frame_index=reference_frame_index)
    print("Writing world positions for",len(star_list),"stars into",init.posdir)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list), total=len(star_list)):
        pass

# TODO check that JD of first line is equal to JD of reference frame !
def world_pos(star, wcs, reference_frame_index):
    #print("star", star)
    f = open(get_pos_filename(star))
    #print("star file opened", star)
    pixel_coords = f.readlines()[2+reference_frame_index].split()[1:3]
    f.close()
    #print("pixel coords read of star", star, pixel_coords)
    world_coords = wcs.all_pix2world(float(pixel_coords[0]), float(pixel_coords[1]), 0, ra_dec_order=True)
    #print("world coords for star", star, world_coords)
    f2 = open(get_worldpos_filename(star), 'w')
    f2.write(str(world_coords[0]) + " " + str(world_coords[1]))
    f2.close()

# helper function
def get_pos_filename(star):
    return init.posdir + "pos_" + str(star).zfill(5) + ".txt"

def get_worldpos_filename(star):
    return init.worldposdir + "worldpos_" + str(star).zfill(5) + ".txt"

def run_determine_reference_frame():
    write_convert_fits()
    write_photometry()

def run_do_rest(reference_phot, do_match, do_munifind, do_lightcurve, do_pos, do_calibrate, do_charts):
    if do_match: write_match(reference_phot)
    if do_munifind:
        write_munifind()
        check_stars_list = do_best_comparison_stars(12)
        with open('check_stars_list.bin', 'wb') as fp:
            pickle.dump(check_stars_list, fp)
        print("check_stars_list: ", check_stars_list)
        write_munifind_check_stars(check_stars_list[0])
    else:
        with open ('check_stars_list.bin', 'rb') as fp:
            check_stars_list = pickle.load(fp)

    if do_lightcurve: do_write_curve(init.star_list, check_stars_list)
    if do_pos: do_write_pos(init.star_list, check_stars_list)
    if do_calibrate:
        wcs = do_calibration.calibrate()
        reference_frame_index = do_calibration.find_reference_in_files(init.fitsdir)
        do_world_pos(wcs, init.star_list, reference_frame_index)
    if do_charts:
        do_charts.run(reference_frame_index)



#logger = mp.log_to_stderr()
#logger.setLevel(mp.SUBDEBUG)
#run_determine_reference_frame()
run_do_rest(init.photometrydir+'phot000001.pht', do_match=False, do_munifind=False, do_lightcurve=False,
            do_pos=True, do_calibrate=True, do_charts=False)

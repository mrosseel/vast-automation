import init
import do_calibration
import do_charts
import do_field_charts
import do_stats_charts
import reading
from reading import trash_and_recreate_dir
from reading import reduce_star_list
import pandas as pd
import multiprocessing as mp
import logging
import tqdm
import os, sys
from functools import partial
from subprocess import call
import subprocess
import pickle
import do_aavso_report
import re


# Munifind fields: INDEX MEAN_MAG STDEV GOODPOINTS
def read_munifind(filename):
    df = pd.read_csv(filename, skiprows=[1], sep=' ')
    df.rename(columns={'INDEX': 'STAR'}, inplace=True)
    print("max goodpoints:", df['GOODPOINTS'].max())
    print("min stdev:", df['STDEV'].min())
    print(df.sort_values('STDEV').head())
    return df


# gets the stars with a maximum of measurements and lowest stdev
def getBestComparisonStars(nrOfStars, filename=init.basedir + 'munifind.txt'):
    print('reading filename', filename)
    df = read_munifind(filename)
    maxPoints = df['GOODPOINTS'].max()
    df = df[df['GOODPOINTS'] > maxPoints * 0.99]
    df_lowest_stdev = df.sort_values('STDEV')
    comparison_stars = df_lowest_stdev.head(nrOfStars)
    #print("Comparison stars:\n", comparison_stars.describe())
    return comparison_stars


def do_best_comparison_stars(nrOfStars):
    bestcomps = getBestComparisonStars(nrOfStars)
    check_stars = []
    for index, row in bestcomps.iterrows():
        # print(row, '\n')
        check_stars.append(int(row['STAR']))
    return check_stars


def join_check_stars(check_stars, exclude_star):
    check_stars = filter(lambda star: star != exclude_star, check_stars)
    check_stars_string = ','.join(map(str, check_stars))
    return check_stars_string

def write_convert_fits():
    print("convert fits files to fts")
    trash_and_recreate_dir(init.convfitsdir)
    os.system('konve ' + init.fitsdir + '*.fit -o ' + init.convfitsdir + 'kout??????.fts')

def write_photometry():
    print("write photometry")
    trash_and_recreate_dir(init.photometrydir)
    os.system('muniphot ' + init.convfitsdir + '*.fts -p muniphot.conf -o ' + init.photometrydir + 'phot??????.pht')

def write_match(to_match_photomotry_file, base_photometry_file):
    # find the number part of to_match_photomotry_file
    m = re.search(r'(\d+)', to_match_photomotry_file)
    numbers = m.group(0)
    # see munimatch.c for arguments of config file
    os.system('munimatch -p ' + init.codedir + 'match.conf ' + base_photometry_file + ' ' + init.photometrydir + to_match_photomotry_file + ' -o ' + init.matchedphotometrydir + 'match'+numbers+'.pht')

def write_munifind(aperture, match_file='match*', quiet=False, canonical=False):
    quiet_switch = '-q ' if quiet else ''
    quiet_null = ' > /dev/null' if quiet else ''
    munifind_file = init.basedir + 'munifind.txt' if canonical else init.basedir + 'munifind_' + str(aperture) + '.txt'
    print("writing munifind for aperture", aperture, "file:", munifind_file, "quiet:", quiet)
    os.system('munifind ' + quiet_switch + '-a ' + str(aperture) + ' ' + munifind_file + ' ' + init.matchedphotometrydir + match_file + quiet_null)
    return munifind_file

def write_munifind_check_stars(check_star, aperture):
    print("write munifind check stars using check star:", check_star)
    os.system('munifind -a ' + str(aperture) + ' ' + ' -c ' + str(check_star) + ' ' + init.basedir + 'munifind.txt ' + init.matchedphotometrydir + 'match*')

def write_lightcurve(star, check_stars_list, aperture):
    check_stars = join_check_stars(check_stars_list, star)
    os.system("munilist -a " + str(aperture) + " -q --object " + str(star) + " -v " + str(
        star) + " -c " + check_stars + " " + init.lightcurvedir + 'curve_' + str(star).zfill(
        5) + ".txt " + init.matchedphotometrydir + 'match*.pht >/dev/null')
    # print("--verbose -a ", str(aperture), " -q --object ", str(star), " -v ", str(star), " -c ", check_stars, (init.lightcurve_dir + 'curve_' + str(star).zfill(5) + ".txt "), (init.basedir+'match*.pht  >/dev/null'))
    # !munilist --verbose -a {str(aperture)} -q --object {str(star)} -v {str(star)} -c {str(8)} {lightcurve_dir + str(star) + ".txt"} {init.basedir+'match*.pht'}

def get_reference_frame(file_limit):
    # Calculate or retrieve reference frame and its index
    try:
        reference_frame, reference_frame_index=reading.read_reference_frame()
    except:
        reference_frame = do_calibration.select_reference_frame(file_limit)
        reference_frame_index = do_calibration.find_file_index(init.convfitsdir, reference_frame, '*.fts')
        lines = [reference_frame, str(reference_frame_index)]
        with open(init.basedir + 'reference_frame.txt', 'w') as f:
            f.write('\n'.join(lines))
    print("Reference frame: {}, index: {}".format(reference_frame, reference_frame_index))
    return [reference_frame, reference_frame_index]

# TODO add check stars to this command?
def write_pos(star, check_stars_list, matched_reference_frame, aperture):
    # start = time.time()
    # check_stars = join_check_stars(check_stars_list, star)
    # os.system("munilist -a " + str(aperture)+ " -q --obj-plot --object "+ str(star)+ " " + get_pos_filename(star) + " " + init.matchedphotometrydir+'match*.pht >/dev/null')
    call("munilist -a " + str(aperture) + " -q --obj-plot --object " + str(star) + " " + reading.get_pos_filename(
        star) + " " + matched_reference_frame + ' >/dev/null', shell=True)
    # end = time.time()


def do_write_pos(star_list, check_stars_list, aperture, is_resume, matched_reference_frame):
    if not is_resume:
        trash_and_recreate_dir(init.posdir)
    else:
        star_list = reduce_star_list(star_list, init.posdir)
    pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
    func = partial(write_pos, check_stars_list=check_stars_list, matched_reference_frame=matched_reference_frame, aperture=aperture)
    print("Writing star positions for", len(star_list), "stars into", init.posdir)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list, 10), total=len(star_list)):
        pass


def do_write_curve(star_list, check_stars_list, aperture, is_resume):
    if not is_resume:
        trash_and_recreate_dir(init.lightcurvedir)
    else:
        star_list = reduce_star_list(star_list, init.lightcurvedir)
    pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
    func = partial(write_lightcurve, check_stars_list=check_stars_list, aperture=aperture)
    print("Writing star lightcurves for", len(star_list), "stars into", init.lightcurvedir)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list, 10), total=len(star_list)):
        pass


def do_world_pos(wcs, star_list, reference_frame_index):
    trash_and_recreate_dir(init.worldposdir)
    print("index", reference_frame_index)
    pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
    func = partial(world_pos, wcs=wcs, reference_frame_index=reference_frame_index)
    print("Writing world positions for", len(star_list), "stars into", init.posdir)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list), total=len(star_list)):
        pass


# TODO check that JD of first line is equal to JD of reference frame !
def world_pos(star, wcs, reference_frame_index):
    f = open(reading.get_pos_filename(star))
    pixel_coords = f.readlines()[2 + reference_frame_index].split()[1:3]
    f.close()
    print("pixel coords read of star", star, pixel_coords)
    world_coords = wcs.all_pix2world(float(pixel_coords[0]), float(pixel_coords[1]), 0, ra_dec_order=True)
    print("world coords for star", star, world_coords)
    f2 = open(reading.get_worldpos_filename(star), 'w')
    f2.write(str(world_coords[0]) + " " + str(world_coords[1]))
    f2.close()

def run_do_rest(do_convert_fits, do_photometry, do_match, do_munifind, do_lightcurve, do_lightcurve_resume, do_pos,
                do_pos_resume,
                do_calibrate, do_ml, do_charting, do_phase_diagram, do_field_charting, do_reporting):

    if do_convert_fits:
        write_convert_fits()

    _, reference_frame_index = get_reference_frame(100)

    if do_photometry:
        write_photometry()

    if do_match:
        pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
        ref_frame = do_calibration.find_reference_photometry(reference_frame_index)
        file_list = reading.get_files_in_dir(init.photometrydir)
        file_list.sort()
        func = partial(write_match, base_photometry_file=ref_frame)
        print("Writing matches for", len(file_list), "stars with reference frame", ref_frame)
        trash_and_recreate_dir(init.matchedphotometrydir)
        for _ in tqdm.tqdm(pool.imap_unordered(func, file_list, 10), total=len(file_list)):
            pass

    if do_munifind:
        aperture = do_calibration.find_optimal_aperture('match??????.pht')
        write_munifind(aperture, quiet=True, canonical=True)
        # we used to do something clever here, but the results are exactly the same as doing the normal thing.
        # a bit too exactly even, but for now we just disable it.
        check_stars_list = do_best_comparison_stars(12)
        with open(init.basedir + 'check_stars_list.bin', 'wb') as fp:
            pickle.dump(check_stars_list, fp)
        # print("check_stars_list: ", check_stars_list)
        # write_munifind_check_stars(check_stars_list[0])
        print("Done writing munifind")
    else:
        with open(init.basedir + 'check_stars_list.bin', 'rb') as fp:
            check_stars_list = pickle.load(fp)
        with open(init.basedir + 'best_aperture.txt', 'r') as fp:
            aperture = round(float(next(fp)),1)
    if do_lightcurve: do_write_curve(init.star_list, check_stars_list, aperture, do_lightcurve_resume)

    if do_pos: do_write_pos(init.star_list, check_stars_list, aperture, do_pos_resume,
                            do_calibration.find_reference_matched(reference_frame_index))

    # is this still needed?? maybe for automatic calibration with astrometry.net
    if do_calibrate:
        wcs = do_calibration.calibrate()
        do_world_pos(wcs, init.star_list, 0)  # pass 0 for reference_frame_index because we only write one position
        df = do_calibration.find_target_star(init.ra_deg, init.dec_deg, 50)
        df.to_csv(init.basedir + 'distances_from_target_star.csv')
        print("calibrate df", df)

    if do_ml:
        import do_upsilon  # do it here because it takes some time at startup
        do_upsilon.run(init.star_list)

    chart_premade = False
    chart_upsilon = False
    chart_vsx = True
    chart_custom = False
    star_descriptions = []
    vsx_star_descriptions = do_calibration.get_vsx_in_field(do_calibration.get_star_descriptions(), 0.01)

    if chart_premade:
        with open(init.basedir + 'star_descriptions_to_chart.bin', 'rb') as fp:
            star_descriptions = pickle.load(fp)
    else:
        if (chart_upsilon):
            # now we generate a list of StarDescriptions, with which all further processing will be done
            star_descriptions = do_calibration.get_candidates(0.1)
        elif(chart_vsx):
            star_descriptions = vsx_star_descriptions
        elif(chart_custom):
            star_descriptions = do_calibration.get_star_descriptions(init.star_list)

        if not chart_vsx: # vsx names are always added, but chart_vsx already has the vsx mappings
            star_descriptions = do_calibration.add_vsx_names_to_star_descriptions(star_descriptions)

        if not chart_upsilon:
            star_descriptions = do_calibration.add_upsilon_data(star_descriptions)

        with open(init.basedir + 'star_descriptions_to_chart.bin', 'wb') as fp:
            pickle.dump(star_descriptions, fp)

    # parse comparison star from munifind.txt
    comparison_star = reading.read_comparison_star()
    print('Read comparison star from munifind.txt: ', comparison_star)

    if 0:
        comp_star_description = do_calibration.get_star_descriptions([comparison_star])
        print('Get comparison star description: ', comp_star_description)
        comparison_stars = do_calibration.add_ucac4_to_star_descriptions(comp_star_description)
        print("Got comparison star description + ucac4 id: ", comparison_stars)
        # add ucac4 to star_descriptions
        star_descriptions_ucac4 = do_calibration.add_ucac4_to_star_descriptions(star_descriptions)

        if do_charting or do_phase_diagram:
            print("starting charting / phase diagrams")
            do_charts.run(star_descriptions_ucac4, comparison_stars, do_charting, do_phase_diagram)

    if do_field_charting:
#        do_field_charts.run_standard_field_charts(vsx_star_descriptions)
        do_stats_charts.plot_cumul_histo_detections()

    #import code
    #code.InteractiveConsole(locals=dict(globals(), **locals())).interact()
    if do_reporting:
        print("AAVSO Reporting with: {} stars".format(len(star_descriptions_ucac4)))
        trash_and_recreate_dir(init.aavsoreportsdir)
        for star in star_descriptions_ucac4:
            do_aavso_report.report(init.aavsoreportsdir, star, comp_star_description[0])


    # logger = mp.log_to_stderr()
    # logger.setLevel(mp.SUBDEBUG)

def interact():
    import code
    code.InteractiveConsole(locals=dict(globals(), **locals())).interact()

if __name__ == '__main__':
    print("Calculating", len(init.star_list), "stars from base dir:", init.basedir, "\nconvert_fits:\t", init.do_convert_fits,
    "\nphotometry:\t", init.do_photometry,
    "\nmatch:\t\t", init.do_match,
    "\nmunifind:\t", init.do_munifind,
    "\nlightcurve:\t", init.do_lightcurve,
    "\nlightcurve_res:\t", init.do_lightcurve_resume,
    "\npos:\t\t", init.do_pos,
    "\npos_resume:\t", init.do_pos_resume,
    "\ncalibrate:\t", init.do_calibrate,
    "\nupsilon:\t", init.do_ml,
    "\ncharting:\t", init.do_charting,
    "\nphasediagram:\t", init.do_phase_diagram,
    "\nfield charts:\t", init.do_field_charts,
    "\nreporting:\t", init.do_reporting)
    print("Press Enter to continue...")
    subprocess.call("read -t 10", shell=True, executable='/bin/bash')
    run_do_rest(init.do_convert_fits, init.do_photometry, init.do_match, init.do_munifind, init.do_lightcurve,
    init.do_lightcurve_resume, init.do_pos, init.do_pos_resume, init.do_calibrate,
    init.do_ml, init.do_charting, init.do_phase_diagram, init.do_field_charts, init.do_reporting)

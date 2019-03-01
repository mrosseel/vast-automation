import init
import do_calibration
import do_charts
import do_charts_field
import do_charts_stats
import do_aperture
import do_lightcurve as dolight
from do_lightcurve import join_check_stars_string
import read_photometry
import reading
from reading import trash_and_recreate_dir
from reading import reduce_star_list
import pandas as pd
import numpy as np
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
import glob
import random


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
    print('getBestComparisonStars: reading ', filename)
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
    for _, row in bestcomps.iterrows():
        # print(row, '\n')
        check_stars.append(int(row['STAR']))
    return check_stars

def write_convert_fits():
    logging.info("Convert fits files to fts")
    trash_and_recreate_dir(init.convfitsdir)
    os.system('konve ' + init.fitsdir + '*.fit -o ' + init.convfitsdir + 'kout??????.fts')

def write_photometry(use_config=True, custom_config=None, custom_wildcard=None, custom_outputfile=None,
                     custom_outputdir=None):
    logging.info(f"Writing photometry, use_config is {use_config}")
    trash_and_recreate_dir(init.photometrydir)
    config_file = ''
    if use_config:
        if custom_config is None:
            config_file = f'-p muniphot.conf'
        else:
            config_file = f'-p {custom_config}'
    config_file = f'-p muniphot.conf' if use_config else ''
    files = init.convfitsdir + '*.fts' if custom_wildcard is None else custom_wildcard
    outputdir = init.photometrydir if custom_outputdir is None else custom_wildcard
    outputfile = 'phot??????.pht' if custom_outputfile is None else custom_outputfile
    os.system(f"muniphot {files} {config_file} -o {outputdir + outputfile}")

def write_match(to_match_photomotry_file, base_photometry_file, use_config=True):
    # find the number part of to_match_photomotry_file
    m = re.search(r'(\d+)', to_match_photomotry_file)
    numbers = m.group(0)
    # see munimatch.c for arguments of config file
    config_file = f'-p {init.codedir + "match.conf"}' if use_config else ''
    os.system(f'munimatch {config_file} {base_photometry_file} {init.photometrydir + to_match_photomotry_file} -o {init.matchedphotometrydir + "match"+numbers+".pht"}')

# percentage is between 0 and 1
def write_munifind_dirfile(match_pattern, percentage=1):
    matched_files = glob.glob(init.matchedphotometrydir+match_pattern)
    desired_length = int(len(matched_files) * percentage)
    np.random.seed(42) # for the same percentage, we always get the same selection
    selected_files = np.random.choice(matched_files, size=desired_length, replace=False).tolist()
    with open(init.aperturedir+'/munifind_filelist.txt', mode='wt', encoding='utf-8') as myfile:
        for lines in selected_files:
            print(lines, file = myfile)
    myfile.close
    return len(selected_files)

def write_munifind(aperture, match_file=False, match_pattern='match*', quiet=False, canonical=False, percentage=1):
    quiet_switch = '-q ' if quiet else ''
    quiet_null = ' > /dev/null' if quiet else ''
    munifind_file = init.basedir + 'munifind.txt' if canonical else init.aperturedir + 'munifind_' + str(aperture) + '.txt'
    selected = None
    if match_file:
        selected = write_munifind_dirfile(match_pattern, percentage)
        sourcefiles = f'-i {init.aperturedir}/munifind_filelist.txt'
    else:
        sourcefiles = init.matchedphotometrydir + match_pattern
    print(f"writing munifind for aperture {aperture}, file: {munifind_file}, quiet: {quiet}, nr_sourcefiles: {selected if selected else 'all'}")
    os.system('munifind ' + quiet_switch + '-a ' + str(aperture) + ' ' + munifind_file + ' ' +  sourcefiles + quiet_null)
    return munifind_file

def write_munifind_check_stars(check_star, aperture):
    print("write munifind check stars using check star:", check_star)
    os.system('munifind -a ' + str(aperture) + ' ' + ' -c ' + str(check_star) + ' ' + init.basedir + 'munifind.txt ' + init.matchedphotometrydir + 'match*')

def write_lightcurve(star, check_stars_list_1, aperture):
    check_stars = join_check_stars_string(check_stars_list_1, star)
    print(check_stars)
    os.system("munilist -a " + str(aperture) + " -q --object " + str(star) + " -v " + str(star) + " -c " + check_stars + " " + init.lightcurvedir + 'curve_' + str(star).zfill(5) + ".txt " + init.matchedphotometrydir + 'match*.pht >/dev/null')
    # print("--verbose -a ", str(aperture), " -q --object ", str(star), " -v ", str(star), " -c ", check_stars, (init.lightcurve_dir + 'curve_' + str(star).zfill(5) + ".txt "), (init.basedir+'match*.pht  >/dev/null'))
    # !munilist --verbose -a {str(aperture)} -q --object {str(star)} -v {str(star)} -c {str(8)} {lightcurve_dir + str(star) + ".txt"} {init.basedir+'match*.pht'}

# TODO add check stars to this command?
def write_pos(star, check_stars_list, matched_reference_frame, aperture):
    # start = time.time()
    # check_stars = join_check_stars(check_stars_list, star)
    # os.system("munilist -a " + str(aperture)+ " -q --obj-plot --object "+ str(star)+ " " + get_pos_filename(star) + " " + init.matchedphotometrydir+'match*.pht >/dev/null')
    call("munilist -a " + str(aperture) + " -q --obj-plot --object " + str(star) + " " + reading.get_pos_filename(star) + " " + matched_reference_frame + ' >/dev/null', shell=True)
    # end = time.time()


def do_write_pos(star_list, check_stars_list, aperture, matched_reference_frame, is_resume):
    if not is_resume:
        trash_and_recreate_dir(init.posdir)
    else:
        star_list = reduce_star_list(star_list, init.posdir)
    pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
    func = partial(write_pos, check_stars_list=check_stars_list, matched_reference_frame=matched_reference_frame, aperture=aperture)
    print("Writing star positions for", len(star_list), "stars into", init.posdir)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list, 10), total=len(star_list)):
        pass

def do_write_curve(star_list, check_stars_1, aperture, is_resume):
    # if not is_resume:
    #     trash_and_recreate_dir(init.lightcurvedir)
    # else:
    #     star_list = reduce_star_list(star_list, init.lightcurvedir)
    #check_stars_0 = np.array(check_stars_list_1) - 1
    pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
    func = partial(write_lightcurve, check_stars_list_1=check_stars_1, aperture=aperture)
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
    pixel_coords = f.readlines()[2].split()[1:3] # there is only one position line, that of the reference frame
    f.close()
    logging.debug(f"pixel coords read of star {star}, {pixel_coords}")
    world_coords = wcs.all_pix2world(float(pixel_coords[0]), float(pixel_coords[1]), 0, ra_dec_order=True)
    logging.debug(f"world coords for star {star}, {world_coords}")
    f2 = open(reading.get_worldpos_filename(star), 'w')
    f2.write(str(world_coords[0]) + " " + str(world_coords[1]))
    f2.close()

def run_do_rest(do_conf_phot, do_conf_match, do_convert_fits, do_photometry, do_match, do_aperture_search, do_lightcurve, do_pos,
                do_ml, do_lightcurve_plot, do_phase_diagram, do_field_charting, do_reporting):
    if do_convert_fits:
        logging.info("Converting fits...")
        write_convert_fits()

    # either read the previous reference frame or calculate a new one
    _, _, reference_frame_index = do_calibration.get_reference_frame(100, do_calibration.select_reference_frame_jpeg)

    # get wcs model from the reference header. Used in writing world positions and field charts
    wcs = do_calibration.get_wcs(init.reference_header)
    comparison_stars_1=-1
    apertures=None
    aperture=None
    apertureidx=None

    if do_photometry:
        logging.info("Writing photometry...")
        write_photometry(use_config=do_conf_phot, custom_wildcard=init.photometry_wildcard)

    if do_match:
        logging.info("Performing matching...")
        pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
        ref_frame = do_calibration.find_reference_photometry(reference_frame_index)
        file_list = reading.get_files_in_dir(init.photometrydir)
        file_list.sort()
        func = partial(write_match, base_photometry_file=ref_frame, use_config=do_conf_match) # Not using config file for testing
        print("Writing matches for", len(file_list), "stars with reference frame", ref_frame)
        trash_and_recreate_dir(init.matchedphotometrydir)
        for _ in tqdm.tqdm(pool.imap_unordered(func, file_list, 10), total=len(file_list)):
            pass

    if do_aperture_search:
        logging.info("Searching best aperture...")
        stddevs, _, apertures, apertureidx, _, _, comparison_stars_1 = do_aperture.main(the_dir=init.matchedphotometrydir, percentage=init.aperture_find_percentage)
        aperture = apertures[apertureidx]
        # with open(init.basedir + 'check_stars_list.bin', 'wb') as fp:
        #     pickle.dump(comparison_stars_1, fp)
        np.savetxt(init.basedir + "comparison_stars_1.txt", comparison_stars_1, fmt='%d', delimiter=';')
        np.savetxt(init.basedir + "apertures.txt", apertures, fmt='%.2f', delimiter=';')
        np.savetxt(init.basedir + "apertureidx_best.txt", [apertureidx], fmt='%d')
        logging.debug("Done writing aperture search results")
    else:
        logging.info("Loading best aperture and compstars...")
        comparison_stars_1, apertures, apertureidx, aperture = reading.aperture_and_compstars()
        logging.info(f"comparison stars: {comparison_stars_1}")
        logging.info(f"aperture: {aperture}, apertures:{apertures}")

    logging.info("Loading photometry...")
    jd, fwhm, star_result = read_photometry.read_photometry(init.star_list, apertureidx)

    # Calculate the quality
    # star_select: jd[nrfiles], fwhm[nrfiles], star_result[nrfiles, nrstars, 2]



    if do_pos:
        logging.info("Writing positions of all stars on the reference image...")
        reference_matched = do_calibration.find_reference_matched(reference_frame_index)
        print("reference match is ", reference_matched)
        do_write_pos(init.star_list, comparison_stars_1, aperture, reference_matched, is_resume=False)
        do_world_pos(wcs, init.star_list, reference_frame_index)

    logging.debug("Before do ml")
    if do_ml:
        logging.debug("Doing ML detection of variable stars...")
        import do_upsilon  # do it here because it takes some time at startup
        do_upsilon.run(init.star_list)

    chart_premade = False
    chart_upsilon = False
    chart_vsx = True
    chart_custom = False
    star_descriptions = []
    logging.info("Getting vsx in field")
    vsx_star_descriptions = do_calibration.get_vsx_in_field(do_calibration.get_star_descriptions(), 0.01)

    if chart_premade:
        logging.info("Loading premade star_descriptions: star_descriptions_to_chart.bin")
        with open(init.basedir + 'star_descriptions_to_chart.bin', 'rb') as fp:
            star_descriptions = pickle.load(fp)
    else:
        if (chart_upsilon):
            # now we generate a list of StarDescriptions, with which all further processing will be done
            logging.info("Setting star_descriptions to upsilon candidates")
            star_descriptions = do_calibration.get_candidates(0.1)
        elif(chart_vsx):
            logging.info("Setting star_descriptions to vsx stars")
            star_descriptions = vsx_star_descriptions
        elif(chart_custom):
            logging.info(f"Setting star_descriptions to custom: {init.star_list}")
            star_descriptions = do_calibration.get_star_descriptions(init.star_list)

        if not chart_vsx: # vsx names are always added, but chart_vsx already has the vsx mappings
            logging.info("Adding vsx names to star_descriptions")
            star_descriptions = do_calibration.add_vsx_names_to_star_descriptions(star_descriptions)

        # if not chart_upsilon:
        #     logging.info("Adding upsilon data to star_descriptions")
        #     star_descriptions = do_calibration.add_upsilon_data(star_descriptions)

        logging.info("Writing star_descriptions_to_chart.bin...")
        with open(init.basedir + 'star_descriptions_to_chart.bin', 'wb') as fp:
            pickle.dump(star_descriptions, fp)

    logging.info("Printing star_descriptions coordinates:")
    # for star in star_descriptions:
    #     print(star.local_id, star.coords, star.match[0].coords, star.match)

    logging.info(f'Getting star description of comparison star: {comparison_stars_1[:1]}')
    comp_star_description = do_calibration.get_star_descriptions(comparison_stars_1[:1])
    logging.info(f'Adding ucac4 info to comparison stars: {comp_star_description}')
    comparison_stars_1_desc = do_calibration.add_ucac4_to_star_descriptions(comp_star_description)

    logging.debug(f"Comparison stars: {comparison_stars_1_desc[0]}")
    if np.isnan(comparison_stars_1_desc[0].vmag):
        print("Comparison star has nan vmag, will screw up everything coming after")
        exit()
    # add ucac4 to star_descriptions (why???)
    logging.info(f"Adding ucac4 to all stars of interest: {star_descriptions}")
    star_descriptions_ucac4 = do_calibration.add_ucac4_to_star_descriptions(star_descriptions)

    if do_lightcurve:
        logging.info(f"Writing lightcurves... {[x.local_id for x in star_descriptions_ucac4]}")
        chosen_stars = [x.local_id for x in star_descriptions_ucac4]
        dolight.write_lightcurves(chosen_stars,
                                  comparison_stars_1, aperture, int(apertureidx), jd, fwhm, star_result)
        # do_write_curve(chosen_stars, comparison_stars_1, aperture, False)

    if do_lightcurve_plot or do_phase_diagram:
        logging.info("starting charting / phase diagrams...")
        do_charts.run(star_descriptions_ucac4, comparison_stars_1_desc, do_lightcurve_plot, do_phase_diagram)

    if do_field_charting:
        logging.info("Starting field chart plotting...")
        do_charts_field.run_standard_field_charts(vsx_star_descriptions, wcs)
        do_charts_stats.main(fwhm)

    #import code
    #code.InteractiveConsole(locals=dict(globals(), **locals())).interact()
    if do_reporting:
        logging.info(f"AAVSO Reporting with: {len(star_descriptions_ucac4)} stars")
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
    "\nUseConf phot:\t", init.do_conf_phot,
    "\nUseConf match:\t", init.do_conf_match,
    "\nphotometry:\t", init.do_photometry,
    "\nmatch:\t\t", init.do_match,
    "\naperture search:\t", init.do_aperture_search,
    "\nlightcurve:\t", init.do_lightcurve,
    "\npos:\t\t", init.do_pos,
    "\nupsilon:\t", init.do_ml,
    "\nlightcurve plot:\t", init.do_lightcurve_plot,
    "\nphasediagram:\t", init.do_phase_diagram,
    "\nfield charts:\t", init.do_field_charts,
    "\nreporting:\t", init.do_reporting)
    print("Press Enter to continue...")
    subprocess.call("read -t 10", shell=True, executable='/bin/bash')
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    run_do_rest(init.do_conf_phot, init.do_conf_match, init.do_convert_fits, init.do_photometry, init.do_match, init.do_aperture_search, init.do_lightcurve,
                init.do_pos, init.do_ml, init.do_lightcurve_plot, init.do_phase_diagram,
                init.do_field_charts, init.do_reporting)

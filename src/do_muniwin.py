import init
import do_calibration
import do_charts
import do_charts_field
import do_charts_stats
import do_aperture
import do_compstars
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

def write_photometry(config_file=init.basedir+'muniphot.conf', files=None, outputfile_prefix='phot',
                         outputdir=init.photometrydir):
    trash_and_recreate_dir(init.photometrydir)
    if files is None:
        files = glob.glob(init.convfitsdir+"*.fts")
    logging.info(f"Writing photometry, config_file:{config_file}, outputdir: {outputdir}, outputfile_prefix: {outputfile_prefix}...")
    config_file_param = f'-p {config_file}' if config_file is not '' else ''
    pool = mp.Pool(init.nr_threads*2, maxtasksperchild=100)
    func = partial(command_caller, config_file_param=config_file_param, outputdir=outputdir, outputfile_prefix=outputfile_prefix)
    print("Writing star lightcurves for", len(files), "stars into", init.lightcurvedir)
    for _ in tqdm.tqdm(pool.imap_unordered(func, files, 50), total=len(files)):
        pass

def command_caller(fits, config_file_param, outputdir, outputfile_prefix):
    if isinstance(fits, list):
        fits = ' '.join(fits)
    command = f"muniphot {fits} {config_file_param} -o {outputdir+outputfile_prefix}{int(''.join(filter(str.isdigit, fits))):06d}.pht"
    logging.debug(f'Photometry command is: {command}')
    subprocess.call(command, shell=True)

def write_match(to_match_photomotry_file, base_photometry_file, to_match_is_full_path=False,
                config_file=init.basedir+'match.conf', inputdir=init.photometrydir, outputdir=init.matchedphotometrydir):
    # find the number part of to_match_photomotry_file
    m = re.search(r'(\d+)(.pht)', to_match_photomotry_file)

    numbers = m.group(0)[:-4]
    print("numbers is",numbers, to_match_photomotry_file, m)
    # see munimatch.c for arguments of config file
    config_file = f'-p {config_file}' if config_file is not '' else ''
    photometry_file_full_path = init.photometrydir + to_match_photomotry_file if to_match_is_full_path is False else to_match_photomotry_file
    command = f'munimatch {config_file} {base_photometry_file} {photometry_file_full_path} -o {outputdir + "match" + numbers + ".pht"}'
    print("Munimatch command:", command)
    subprocess.call(command, shell=True)

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

def do_write_curve_old(star_list, check_stars_1, aperture, is_resume):
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

# experimental multithreading of new lightcurve writing method, not yet tested/used
def do_write_curve(star_list, check_stars_1, aperture, apertureidx, jd, fwhm, star_result_):
    pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
    func = partial(dolight.write_lightcurves, check_stars_list_1=check_stars_1, aperture=aperture,
                   apertureidx=apertureidx, jd=jd, fwhm=fwhm, star_result_=star_result_)
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

def run_do_rest(do_convert_fits, do_photometry, do_match, do_aperture_search, do_lightcurve, do_pos,
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
        logging.info("Writing photometry with config file {init.conf_phot}...")
        write_photometry(config_file=init.conf_phot)

    if do_match:
        logging.info("Performing matching...")
        pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
        ref_frame = do_calibration.find_reference_photometry(reference_frame_index)
        file_list = reading.get_files_in_dir(init.photometrydir)
        file_list.sort()
        func = partial(write_match, base_photometry_file=ref_frame)
        print("Writing matches for", len(file_list), "stars with reference frame", ref_frame)
        trash_and_recreate_dir(init.matchedphotometrydir)
        for _ in tqdm.tqdm(pool.imap_unordered(func, file_list, 10), total=len(file_list)):
            pass

    if do_aperture_search:
        logging.info("Searching best aperture...")
        # getting aperture
        stddevs, _, apertures, apertureidx, _, _, counts = do_aperture.main(the_dir=init.matchedphotometrydir, percentage=init.aperture_find_percentage)
        aperture = apertures[apertureidx]
        # getting compstars
        comparison_stars_1, comparison_stars_1_desc = do_compstars.select_compstars(stddevs, apertureidx, counts)
        logging.info(f"Comparison stars_1: {comparison_stars_1}")
        # saving all calculated data
        np.savetxt(init.basedir + "comparison_stars_1.txt", comparison_stars_1, fmt='%d', delimiter=';')
        np.savetxt(init.basedir + "apertures.txt", apertures, fmt='%.2f', delimiter=';')
        np.savetxt(init.basedir + "apertureidx_best.txt", [apertureidx], fmt='%d')
        with open(init.basedir + 'comparison_stars_1_desc.bin', 'wb') as compfile:
            pickle.dump(comparison_stars_1_desc, compfile)
        logging.debug("Done writing aperture search results")
    else:
        logging.info("Loading best aperture and compstars...")
        comparison_stars_1, comparison_stars_1_desc, apertures, apertureidx, aperture = reading.read_aperture_and_compstars()
        logging.info(f"comparison stars: {comparison_stars_1}")
        logging.info(f"aperture: {aperture}, apertures:{apertures}")

    logging.info("Loading photometry...")
    jd, fwhm, nrstars, star_result = read_photometry.read_photometry(init.star_list, apertureidx)

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

    # add ucac4 to star_descriptions for AAVSO reporting (specifically for unknown stars)
    # logging.info(f"Adding ucac4 to all stars of interest: {star_descriptions}")
    star_descriptions_ucac4 = do_calibration.add_ucac4_to_star_descriptions(star_descriptions)

    if do_lightcurve:
        logging.info(f"Writing lightcurves... {[x.local_id for x in star_descriptions_ucac4]}")
        chosen_stars = [x.local_id for x in star_descriptions_ucac4]
        dolight.write_lightcurves(chosen_stars,
                                  comparison_stars_1, aperture, int(apertureidx), jd, fwhm, star_result)
        # do_write_curve(chosen_stars,
        #                           comparison_stars_1, aperture, int(apertureidx), jd, fwhm, star_result)
        #
        # do_write_curve(chosen_stars, comparison_stars_1, aperture, False)

    if do_lightcurve_plot or do_phase_diagram:
        logging.info("starting charting / phase diagrams...")
        print("comparison stars decs:", comparison_stars_1_desc)
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
            do_aavso_report.report(init.aavsoreportsdir, star, comparison_stars_1_desc[0])


    # logger = mp.log_to_stderr()
    # logger.setLevel(mp.SUBDEBUG)

def interact():
    import code
    code.InteractiveConsole(locals=dict(globals(), **locals())).interact()

if __name__ == '__main__':
    print("Calculating", len(init.star_list), "stars from base dir:", init.basedir, "\nconvert_fits:\t", init.do_convert_fits,
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
    run_do_rest(init.do_convert_fits, init.do_photometry, init.do_match, init.do_aperture_search, init.do_lightcurve,
                init.do_pos, init.do_ml, init.do_lightcurve_plot, init.do_phase_diagram,
                init.do_field_charts, init.do_reporting)

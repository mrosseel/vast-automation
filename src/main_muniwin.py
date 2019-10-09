import do_calibration
from importlib import reload
import do_lightcurve
do_lightcurve = reload(do_lightcurve)
import do_charts
import do_charts_field
import do_aperture
import do_compstars
import read_photometry
import reading
reading = reload(reading)
from reading import trash_and_recreate_dir
from reading import reduce_star_list
from star_description import StarDescription
from photometry_blob import PhotometryBlob
import numpy as np
import multiprocessing as mp
import tqdm
import os
from functools import partial
from subprocess import call
import subprocess
import pickle
import do_aavso_report
import re
import glob
import utils
import logging
from init_loader import init, settings

def write_convert_fits():
    logging.info("Convert fits files to fts")
    trash_and_recreate_dir(settings.convfitsdir)
    os.system('konve ' + settings.fitsdir + '*.fit -o ' + settings.convfitsdir + 'kout??????.fts')


def write_photometry(config_file=settings.basedir + 'muniphot.conf', files=None, outputfile_prefix='phot',
                     outputdir=settings.photometrydir):
    trash_and_recreate_dir(settings.photometrydir)
    if files is None:
        files = glob.glob(settings.convfitsdir + "*.fts")
    logging.info(
        f"Writing photometry, config_file:{config_file}, outputdir: {outputdir}, outputfile_prefix: {outputfile_prefix}...")
    config_file_param = f'-p {config_file}' if config_file is not '' else ''
    pool = mp.Pool(init.nr_threads * 2, maxtasksperchild=100)
    func = partial(command_caller, config_file_param=config_file_param, outputdir=outputdir,
                   outputfile_prefix=outputfile_prefix)
    logging.info(f"Writing star photometry for {len(files)} stars into {outputdir}")
    logging.info(files)
    for _ in tqdm.tqdm(pool.imap_unordered(func, files, 5), total=len(files)):
        pass


def write_match(to_match_photomotry_file, base_photometry_file, to_match_is_full_path=False,
                config_file=settings.conf_match, inputdir=settings.photometrydir,
                outputdir=settings.matchedphotometrydir):
    # find the number part of to_match_photomotry_file
    m = re.search(r'(\d+)(.pht)', to_match_photomotry_file)

    numbers = m.group(0)[:-4]
    logging.info(f"numbers is {numbers} {to_match_photomotry_file} {m}")
    # see munimatch.c for arguments of config file
    config_file = f'-p {config_file}' if config_file is not '' else ''
    photometry_file_full_path = settings.photometrydir + to_match_photomotry_file if to_match_is_full_path is False else to_match_photomotry_file
    command = f'munimatch {config_file} {base_photometry_file} {photometry_file_full_path} -o {outputdir + "match" + numbers + ".pht"}'
    logging.debug(f"Munimatch command: {command}")
    subprocess.call(command, shell=True)


# TODO read the x y position from the photometry files, ditch munilist
def write_pos(star, matched_reference_frame, aperture):
    call("munilist -a " + str(aperture) + " -q --obj-plot --object " + str(star) + " "
         + get_pos_filename(star) + " " + matched_reference_frame + ' >/dev/null', shell=True)


def do_write_pos(star_list, aperture, matched_reference_frame, is_resume):
    if not is_resume:
        trash_and_recreate_dir(settings.posdir)
    else:
        star_list = reduce_star_list(star_list, settings.posdir)
    pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
    func = partial(write_pos, matched_reference_frame=matched_reference_frame,
                   aperture=aperture)
    logging.info(f"Writing star positions for {len(star_list)} stars into {settings.posdir}")
    logging.info(dir(init))
    logging.info(dir(settings))
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list, 10), total=len(star_list)):
        pass


def do_world_pos(wcs, star_list, reference_frame_index):
    trash_and_recreate_dir(settings.worldposdir)
    logging.info(f"index {reference_frame_index}")
    pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
    func = partial(world_pos, wcs=wcs, reference_frame_index=reference_frame_index)
    logging.info(f"Writing world positions for {len(star_list)} stars into {settings.posdir}")
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list), total=len(star_list)):
        pass


# TODO check that JD of first line is equal to JD of reference frame !
def world_pos(star, wcs, reference_frame_index):
    f = open(get_pos_filename(star))
    pixel_coords = f.readlines()[2].split()[1:3]  # there is only one position line, that of the reference frame
    f.close()
    logging.debug(f"pixel coords read of star {star}, {pixel_coords}")
    world_coords = wcs.all_pix2world(float(pixel_coords[0]), float(pixel_coords[1]), 0, ra_dec_order=True)
    logging.debug(f"world coords for star {star}, {world_coords}")
    f2 = open(get_worldpos_filename(star), 'w')
    f2.write(str(world_coords[0]) + " " + str(world_coords[1]))
    f2.close()

# helper function
def get_pos_filename(star):
    return settings.posdir + "pos_" + str(star).zfill(5) + ".txt"

def get_worldpos_filename(star):
    return settings.worldposdir + "worldpos_" + str(star).zfill(5) + ".txt"


def command_caller(fits, config_file_param, outputdir, outputfile_prefix):
    if isinstance(fits, list):
        fits = ' '.join(fits)
    logging.info(f"'{fits}' '{config_file_param}'")
    command = f"muniphot {fits} {config_file_param} -o {outputdir + outputfile_prefix}{int(''.join(filter(str.isdigit, fits))):06d}.pht"
    logging.debug(f'Photometry command is: {command}')
    subprocess.call(command, shell=True)


def run_do_rest(do_convert_fits, do_photometry, do_match, do_compstars_flag, do_aperture_search, do_lightcurve_flag, do_pos,
                do_ml, do_lightcurve_plot, do_phase_diagram, do_field_charting, do_reporting, args):
    if do_convert_fits:
        logging.info("Converting fits...")
        write_convert_fits()

    # either read the previous reference frame or calculate a new one
    _, _, reference_frame_index = do_calibration.get_reference_frame(100, do_calibration.select_reference_frame_jpeg)

    logging.info(f"reference header is {settings.reference_header}")
    # get wcs model from the reference header. Used in writing world positions and field charts
    wcs = do_calibration.get_wcs(settings.reference_header)
    apertures = None
    aperture = None
    apertureidx = None
    comparison_stars_1, comparison_stars_1_desc = None, None
    photometry_blob = PhotometryBlob() # not yet used

    if do_photometry:
        logging.info(f"Writing photometry with config file {settings.conf_phot}...")
        write_photometry(config_file=settings.conf_phot)

    if do_match:
        logging.info("Performing matching...")
        pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
        ref_frame = do_calibration.find_reference_photometry(reference_frame_index)
        file_list = utils.get_files_in_dir(settings.photometrydir)
        file_list.sort()
        func = partial(write_match, base_photometry_file=ref_frame)
        logging.info(f"Writing matches for {len(file_list)} stars with reference frame {ref_frame}")
        trash_and_recreate_dir(settings.matchedphotometrydir)
        for _ in tqdm.tqdm(pool.imap_unordered(func, file_list, 10), total=len(file_list)):
            pass

    # if we need to calculate aperture or comparison stars, we need the photometry blob
    if do_aperture_search or init.comparison_stars is None:
        photometry_blob = do_aperture.get_photometry_blob(the_dir=settings.matchedphotometrydir,
                                                          percentage=init.aperture_find_percentage)

    if do_aperture_search:
        logging.info("Searching best aperture...")
        # getting aperture
        stddevs = None
        counts = None
        apertures = [x for x in do_aperture.get_apertures()]

        if init.aperture is None:
            logging.error("automatic aperture search is not supported yet !!!")
        else:
            apertureidx = np.abs(np.array(apertures) - init.aperture).argmin()

        aperture = apertures[apertureidx]
        # saving all calculated data
        np.savetxt(settings.basedir + "apertures.txt", apertures, fmt='%.2f', delimiter=';')
        np.savetxt(settings.basedir + "apertureidx_best.txt", [apertureidx], fmt='%d')
        logging.debug("Done writing aperture search results")
    else:
        logging.info("Loading best aperture...")
        apertures, apertureidx, aperture = reading.read_aperture()
        logging.info(f"aperture: {aperture}, apertures:{apertures}")

    if do_pos:
        logging.info("Writing positions of all stars on the reference image...")
        reference_matched = do_calibration.find_reference_matched(reference_frame_index)
        logging.info(f"reference match is {reference_matched}")
        do_write_pos(init.star_list, aperture, reference_matched, is_resume=False)
        do_world_pos(wcs, init.star_list, reference_frame_index)

    # get comparison stars, needs world positions !!!
    if init.comparison_stars is None:
        comparison_stars_1, comparison_stars_1_desc = do_compstars.get_calculated_compstars(apertureidx, photometry_blob)
    else:
        comparison_stars_1, comparison_stars_1_desc = do_compstars.get_fixed_compstars()

    if do_ml:
        logging.info("Doing ML detection of variable stars...")
        import do_upsilon  # do it here because it takes some time at startup
        do_upsilon.run(init.star_list)

    if (do_lightcurve_flag or do_field_charting):
        logging.info("Loading photometry...")
        jd, fwhm, nrstars, star_result = read_photometry.read_photometry(init.star_list, apertureidx)


    # costruction of the star descriptions list
    star_descriptions = construct_star_descriptions(args, do_compstars_flag, comparison_stars_1, comparison_stars_1_desc)

    selected_filter = partial(utils.catalog_filter, catalog_name='SELECTED')
    selected_stars = list(filter(selected_filter, star_descriptions))

    logging.debug(f"Selecting chosen stars with {len(selected_stars)} stars selected.")

    if do_lightcurve_flag:
        logging.info(f"Writing lightcurves...")
        do_lightcurve.write_lightcurves([x.local_id for x in star_descriptions],
                                        comparison_stars_1, aperture, int(apertureidx), jd, fwhm, star_result)

    if do_lightcurve_plot or do_phase_diagram:
        logging.info("starting charting / phase diagrams...")
        do_charts.run(star_descriptions, comparison_stars_1_desc, do_lightcurve_plot, do_phase_diagram)

    if do_field_charting:
        logging.info("Starting field chart plotting...")
        do_charts_field.run_standard_field_charts(selected_stars, wcs, settings.fieldchartsdirs, settings.reference_header)

    # import code
    # code.InteractiveConsole(locals=dict(globals(), **locals())).interact()
    if do_reporting:
        # star_descriptions_ucac4 = do_calibration.add_ucac4_to_star_descriptions(star_descriptions)
        logging.info(f"AAVSO Reporting with: {len(selected_stars)} stars")
        trash_and_recreate_dir(settings.aavsoreportsdir)
        for star in selected_stars:
            do_aavso_report.report(star, settings.aavsoreportsdir, comparison_stars_1_desc[0], filter=None)

# make a list of star_descriptions containing all selected stars
def construct_star_descriptions(args, do_compstars_flag, comparison_stars_1, comparison_stars_1_desc):
    # Start with the list of all detected stars
    star_descriptions = do_calibration.get_star_descriptions(init.star_list)
    if args.laststars:
        file_to_load = settings.basedir + 'star_descriptions_to_chart.bin'
        logging.info(f"Loading premade star_descriptions: {file_to_load}")
        with open(file_to_load, 'rb') as fp:
            star_descriptions = pickle.load(fp)
    else:
        if args.upsilon:
            # now we generate a list of StarDescriptions, with which all further processing will be done
            logging.info("Setting star_descriptions to upsilon candidates")
            star_descriptions = do_calibration.add_candidates_to_star_descriptions(star_descriptions, 0.1)

        star_descriptions, results_ids_0 = do_calibration.add_vsx_names_to_star_descriptions(star_descriptions, settings.vsxcatalogdir, 0.01)
        results_ids_0.sort()

        # write the vsx stars used into a file
        with open(settings.basedir + 'vsx_stars.txt', 'wt') as fp:
            for vsx_id in results_ids_0:
                fp.write(f"{vsx_id+1}:\t{star_descriptions[vsx_id].aavso_id}\t{utils.get_hms_dms(star_descriptions[vsx_id].coords)}\n")

        if args.vsx:
            do_calibration.log_star_descriptions(star_descriptions, results_ids_0)
            do_calibration.add_selected_match_to_stars(star_descriptions, results_ids_0)  # select star ids

        if args.starfile:
            with open(settings.basedir + args.starfile, 'r') as fp:
                lines = fp.readlines()
                starlist_1 = [x.rstrip() for x in lines]
                logging.debug(f"The list of stars read from the starfile is: {starlist_1} ")
                starlist_1 = [int(x) for x in filter(str.isdigit, starlist_1)]
                logging.debug(f"The list of stars read from the starfile is: {starlist_1} ")
                logging.info(f"Selecting {starlist_1} stars added by {args.starfile}")
                starlist_0 = [x-1 for x in starlist_1]
                do_calibration.add_selected_match_to_stars(star_descriptions, starlist_0)  # select star ids

        # compstar data is added to star descriptions ==> no it's not
        if do_compstars_flag:
            logging.info("Getting comparison stars...")
            logging.info(f"Comparison stars_1: {comparison_stars_1}")
            np.savetxt(settings.basedir + "comparison_stars_1.txt", comparison_stars_1, fmt='%d', delimiter=';')
            with open(settings.basedir + 'comparison_stars_1_desc.bin', 'wb') as compfile:
                pickle.dump(comparison_stars_1_desc, compfile)
        else:
            logging.info("Loading best compstars...")
            comparison_stars_1, comparison_stars_1_desc = reading.read_compstars()
            logging.info(f"comparison stars: {comparison_stars_1}")

        logging.info("Writing star_descriptions_to_chart.bin...")
        with open(settings.basedir + 'star_descriptions_to_chart.bin', 'wb') as fp:
            pickle.dump(star_descriptions, fp)
    return star_descriptions


def interact():
    import code
    code.InteractiveConsole(locals=dict(globals(), **locals())).interact()


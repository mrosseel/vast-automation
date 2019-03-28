import init
import do_calibration
import do_charts
import do_charts_field
import do_charts_stats
import do_aperture
import do_compstars
import do_lightcurve as dolight
import read_photometry
import reading
from reading import trash_and_recreate_dir
from reading import reduce_star_list
import numpy as np
import multiprocessing as mp
import logging
import tqdm
import os
from functools import partial
from subprocess import call
import subprocess
import pickle
import do_aavso_report
import re
import glob
import argparse

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
    logging.debug(f"Munimatch command: {command}")
    subprocess.call(command, shell=True)


# TODO read the x y position from the photometry files, ditch munilist
def write_pos(star, check_stars_list, matched_reference_frame, aperture):
    call("munilist -a " + str(aperture) + " -q --obj-plot --object " + str(star) + " " + reading.get_pos_filename(star) + " " + matched_reference_frame + ' >/dev/null', shell=True)

def do_write_pos(star_list, check_stars_list, aperture, matched_reference_frame, is_resume):
    if not is_resume:
        trash_and_recreate_dir(init.posdir)
    else:
        star_list = reduce_star_list(star_list, init.posdir)
    pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
    func = partial(write_pos, check_stars_list=check_stars_list, matched_reference_frame=matched_reference_frame, aperture=aperture)
    logging.info(f"Writing star positions for {len(star_list)} stars into {init.posdir}")
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

def command_caller(fits, config_file_param, outputdir, outputfile_prefix):
    if isinstance(fits, list):
        fits = ' '.join(fits)
    command = f"muniphot {fits} {config_file_param} -o {outputdir+outputfile_prefix}{int(''.join(filter(str.isdigit, fits))):06d}.pht"
    logging.debug(f'Photometry command is: {command}')
    subprocess.call(command, shell=True)


def run_do_rest(do_convert_fits, do_photometry, do_match, do_compstars_flag, do_aperture_search, do_lightcurve, do_pos,
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
        stddevs = None
        counts = None
        # stddevs, _, apertures, apertureidx, _, _, counts = do_aperture.main(the_dir=init.matchedphotometrydir, percentage=init.aperture_find_percentage)
        apertures = [x for x in do_aperture.get_apertures()]
        apertureidx = np.abs(np.array(apertures) - init.aperture).argmin()
        aperture = apertures[apertureidx]
        # saving all calculated data
        np.savetxt(init.basedir + "apertures.txt", apertures, fmt='%.2f', delimiter=';')
        np.savetxt(init.basedir + "apertureidx_best.txt", [apertureidx], fmt='%d')
        logging.debug("Done writing aperture search results")
    else:
        logging.info("Loading best aperture...")
        apertures, apertureidx, aperture = reading.read_aperture()
        logging.info(f"aperture: {aperture}, apertures:{apertures}")

    if do_pos:
        logging.info("Writing positions of all stars on the reference image...")
        reference_matched = do_calibration.find_reference_matched(reference_frame_index)
        print("reference match is ", reference_matched)
        do_write_pos(init.star_list, comparison_stars_1, aperture, reference_matched, is_resume=False)
        do_world_pos(wcs, init.star_list, reference_frame_index)

    if do_compstars_flag:
        logging.info("Getting comparison stars...")
        # getting compstars
        comparison_stars_1, comparison_stars_1_desc = do_compstars.get_fixed_compstars()
        logging.info(f"Comparison stars_1: {comparison_stars_1}")
        np.savetxt(init.basedir + "comparison_stars_1.txt", comparison_stars_1, fmt='%d', delimiter=';')
        with open(init.basedir + 'comparison_stars_1_desc.bin', 'wb') as compfile:
            pickle.dump(comparison_stars_1_desc, compfile)
    else:
        logging.info("Loading best compstars...")
        comparison_stars_1, comparison_stars_1_desc = reading.read_compstars()
        logging.info(f"comparison stars: {comparison_stars_1}")

    if(do_lightcurve or do_field_charting):
        logging.info("Loading photometry...")
        jd, fwhm, nrstars, star_result = read_photometry.read_photometry(init.star_list, apertureidx)

    # Calculate the quality
    # star_select: jd[nrfiles], fwhm[nrfiles], star_result[nrfiles, nrstars, 2]


    ## Star description construction

    chart_premade = False
    chart_only_upsilon = False
    chart_only_vsx = False
    chart_only_custom = True
    star_descriptions = []

    if chart_premade:
        logging.info("Loading premade star_descriptions: star_descriptions_to_chart.bin")
        with open(init.basedir + 'star_descriptions_to_chart.bin', 'rb') as fp:
            star_descriptions = pickle.load(fp)
    else:
        if (chart_only_upsilon):
            # now we generate a list of StarDescriptions, with which all further processing will be done
            logging.info("Setting star_descriptions to upsilon candidates")
            star_descriptions = do_calibration.get_candidates(0.1)
        elif(chart_only_vsx):
            logging.info("Setting star_descriptions to vsx stars")
            star_descriptions = do_calibration.get_vsx_in_field(do_calibration.get_star_descriptions(), 0.01)
        elif(chart_only_custom):
            logging.info(f"Setting star_descriptions to custom: {init.star_list}")
            star_descriptions = do_calibration.get_star_descriptions(init.star_list)

        if not chart_only_vsx: # vsx names are always added, but chart_only_vsx already has the vsx mappings
            logging.info("Adding vsx names to star_descriptions")
            star_descriptions = do_calibration.add_vsx_names_to_star_descriptions(star_descriptions)


        logging.info("Writing star_descriptions_to_chart.bin...")
        with open(init.basedir + 'star_descriptions_to_chart.bin', 'wb') as fp:
            pickle.dump(star_descriptions, fp)

    # add ucac4 to star_descriptions for AAVSO reporting (specifically for unknown stars)
    # logging.info(f"Adding ucac4 to all stars of interest: {star_descriptions}")
    #star_descriptions_ucac4 = do_calibration.add_ucac4_to_star_descriptions(star_descriptions)

    if do_lightcurve:
        logging.info(f"Writing lightcurves...")
        chosen_stars = [x.local_id for x in star_descriptions]
        dolight.write_lightcurves(chosen_stars,
                                  comparison_stars_1, aperture, int(apertureidx), jd, fwhm, star_result)

    if do_ml:
        logging.info("Doing ML detection of variable stars...")
        import do_upsilon  # do it here because it takes some time at startup
        do_upsilon.run(init.star_list)
    else:
        if not chart_only_upsilon:
            logging.info("Adding upsilon data to star_descriptions")
            star_descriptions = do_calibration.add_upsilon_data(star_descriptions)

    if do_lightcurve_plot or do_phase_diagram:
        logging.info("starting charting / phase diagrams...")
        do_charts.run(star_descriptions, comparison_stars_1_desc, do_lightcurve_plot, do_phase_diagram)

    if do_field_charting:
        logging.info("Starting field chart plotting...")
        do_charts_field.run_standard_field_charts(vsx_star_descriptions, wcs)
        do_charts_stats.main(fwhm)

    #import code
    #code.InteractiveConsole(locals=dict(globals(), **locals())).interact()
    if do_reporting:
        star_descriptions_ucac4 = do_calibration.add_ucac4_to_star_descriptions(star_descriptions)
        logging.info(f"AAVSO Reporting with: {len(star_descriptions_ucac4)} stars")
        trash_and_recreate_dir(init.aavsoreportsdir)
        for star in star_descriptions_ucac4:
            do_aavso_report.report(init.aavsoreportsdir, star, comparison_stars_1_desc[0])



def interact():
    import code
    code.InteractiveConsole(locals=dict(globals(), **locals())).interact()

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    # logger = mp.log_to_stderr()
    # logger.setLevel(mp.SUBDEBUG)

    parser = argparse.ArgumentParser(description='munipack automation cli')
    parser.add_argument('-c', '--chart', help="Generate lightcurve, lightcurve plot and phase diagram plot", nargs='+')
    parser.add_argument('-n', '--nowait', help="Don't wait 10 secs before starting", action="store_true")
    args = parser.parse_args()
    if args.chart:
        print("in chart part")
        logging.info(f"Writing lightcurves... {[x.local_id for x in star_descriptions_ucac4]}")
        chosen_stars = [x.local_id for x in star_descriptions_ucac4]
        dolight.write_lightcurves(chosen_stars,
                                  comparison_stars_1, aperture, int(apertureidx), jd, fwhm, star_result)
        logging.info("starting charting / phase diagrams...")
        print("comparison stars decs:", comparison_stars_1_desc)
        do_charts.run(star_descriptions_ucac4, comparison_stars_1_desc, do_lightcurve_plot, do_phase_diagram)

    else:
        print("Calculating", len(init.star_list), "stars from base dir:", init.basedir, "\nconvert_fits:\t", init.do_convert_fits,
              "\nphotometry:\t", init.do_photometry,
              "\nmatch:\t\t", init.do_match,
              "\naperture:\t", init.do_aperture_search,
              "\npos:\t\t", init.do_pos,
              "\ncompstars:\t", init.do_compstars,
              "\nlightcurve:\t", init.do_lightcurve,
              "\nupsilon:\t", init.do_ml,
              "\nlight plot:\t", init.do_lightcurve_plot,
              "\nphasediagram:\t", init.do_phase_diagram,
              "\nfield charts:\t", init.do_field_charts,
              "\nreporting:\t", init.do_reporting)
        if not args.nowait:
            print("Press Enter to continue...")
            subprocess.call("read -t 10", shell=True, executable='/bin/bash')
        run_do_rest(init.do_convert_fits, init.do_photometry, init.do_match, init.do_compstars, init.do_aperture_search, init.do_lightcurve,
                    init.do_pos, init.do_ml, init.do_lightcurve_plot, init.do_phase_diagram,
                    init.do_field_charts, init.do_reporting)

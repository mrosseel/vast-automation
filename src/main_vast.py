from multiprocessing import cpu_count
import do_calibration
from importlib import reload
import do_charts_vast
import reading
reading = reload(reading)
from star_description import StarDescription
from photometry_blob import PhotometryBlob
import logging
from init_loader import init, settings
from reading import file_selector
import utils
from astropy.coordinates import SkyCoord
from typing import List
import re
import os.path
from functools import partial

vsx_catalog_name = "vsx_catalog.bin"
vsxcatalogdir = '/home/jovyan/work/' + vsx_catalog_name


def run_do_rest(args):

    vastdir = utils.add_trailing_slash(args.datadir)
    logging.info(f"Directory where all files will be read & written: {vastdir}")
    wcs_file = vastdir+'new-image.fits'
    reference_frame = extract_reference_frame(vastdir)
    logging.info(f"The reference frame is {reference_frame}")
    logging.info(f"reference header is {wcs_file}")
    while not os.path.isfile(wcs_file):
        logging.info("Please provide the reference header '{wcs_file}', which is an astrometry.net plate-solve of {reference_frame} and press Enter to continue...")
        subprocess.call("read -t 10", shell=True, executable='/bin/bash')

    selected_files = file_selector(the_dir=vastdir, match_pattern="*.dat")
    nr_selected = len(selected_files)

    logging.info(f"Selected {nr_selected} light curves from vast dir.")

    # get wcs model from the reference header. Used in writing world positions and field charts
    wcs = do_calibration.get_wcs(wcs_file)
    all_stardict = read_stardict(vastdir)
    #print(sorted(all_stardict.keys())[:10])
    logging.info(f"Number of found lightcurves: {len(selected_files)}")
    star_descriptions = construct_star_descriptions(vastdir, wcs, all_stardict, selected_files, args)
    write_augmented_autocandidates(vastdir, star_descriptions)
    write_augmented_all_stars(vastdir, star_descriptions)
    if args.phaseall:
        do_charts_vast.run(star_descriptions, vastdir+'phase/', vastdir+'chart/', cpu_count())
    elif args.candidates:
        candidate_ids = get_autocandidates(vastdir)
        logging.info(f"Selecting {len(candidate_ids)} candidates to plot")
        selected_stars = do_calibration.select_star_descriptions(candidate_ids, star_descriptions)
        do_charts_vast.run(selected_stars, vastdir+'phase/', vastdir+'chart/', cpu_count())


    comparison_stars_1, comparison_stars_1_desc = None, None
    photometry_blob = PhotometryBlob() # not yet used

    # select a subset of star_descriptions, as specified in starfile
    # if args.starfile:
    #     selected_filter = partial(catalog_filter, catalog_name='SELECTED')
    #     selected_stars = list(filter(selected_filter, star_descriptions))
    #     logging.debug(f"Light curve: selecting chosen stars with {len(selected_stars)} stars selected.")
    # else:
    #     logging.debug(f"No starfile, selecting all stars")
    #     selected_stars = star_descriptions
    #

    # if do_lightcurve_plot or do_phase_diagram:
    #     logging.info("starting charting / phase diagrams...")
    #     do_charts_vast.run(star_descriptions, do_lightcurve_plot, do_phase_diagram)

    # if do_field_charting:
    #     logging.info("Starting field chart plotting...")
    #     do_charts_field.run_standard_field_charts(selected_stars, wcs)

    # import code
    # code.InteractiveConsole(locals=dict(globals(), **locals())).interact()
    # if do_reporting:
    #     # star_descriptions_ucac4 = do_calibration.add_ucac4_to_star_descriptions(star_descriptions)
    #     logging.info(f"AAVSO Reporting with: {len(selected_stars)} stars")
    #     trash_and_recreate_dir(settings.aavsoreportsdir)
    #     for star in selected_stars:
    #         do_aavso_report.report(settings.aavsoreportsdir, star, comparison_stars_1_desc[0])


def extract_reference_frame(from_dir):
    # Ref.  image: 2458586.50154 13.04.2019 00:00:41   ../../inputfiles/TXCar/fits/TXCar#45V_000601040_FLAT.fit
    result = [re.findall(r'Ref.  image: (?:.*) (.*)',line) for line in open(from_dir+'vast_summary.log')]
    return [x for x in result if x != []][0][0]

def read_stardict(vastdir):
    stardict = {}
    for line in open(vastdir+'vast_list_of_all_stars.log'):
        splitline = line.split()
        stardict[int(splitline[0])] = (splitline[1], splitline[2])
    return stardict


def get_starid_from_outfile(outfile):
    m = re.search('out(.*).dat', outfile)
    return int(m.group(1).lstrip('0'))

def write_augmented_autocandidates(dir, stars: List[StarDescription]):
    origname = 'vast_autocandidates.log'
    newname = 'vast_autocandidates_pos.txt'
    logging.info(f"Writing {newname}...")
    cachedict = {}
    for sd in stars:
        cachedict[sd.local_id] = sd
    with open(dir+origname, 'r', encoding='utf-8') as infile, open(dir+newname, 'w') as outfile:
        for line in infile:
            linetext = line.rstrip()
            star_id = get_starid_from_outfile(linetext)
            cacheentry = cachedict[star_id]
            outfile.write(f"{linetext}{'' if cacheentry.path is not '' else '*'}\t{cacheentry.aavso_id}\t{utils.get_hms_dms(cacheentry.coords)}\n")


def get_autocandidates(dir) -> List[int]:
    origname = 'vast_autocandidates.log'
    result = []
    with open(dir+origname, 'r', encoding='utf-8') as infile:
        for line in infile:
            linetext = line.rstrip()
            star_id = get_starid_from_outfile(linetext)
            result.append(star_id)
    return result

def write_augmented_all_stars(dir, stars: List[StarDescription]):
    origname = 'vast_list_of_all_stars.log'
    newname = 'vast_list_of_all_stars_pos.txt'
    logging.info(f"Writing {newname}...")
    cachedict = {}
    for sd in stars:
        cachedict[sd.local_id] = sd
    with open(dir+origname, 'r', encoding='utf-8') as infile, open(dir+newname, 'w') as outfile:
        for line in infile:
            star_id = line.split()[0]
            cacheentry = cachedict[int(star_id)]
            outfile.write(f"{star_id}\t{cacheentry.aavso_id}\t{utils.get_hms_dms(cacheentry.coords)}\n")

def write_vsx_stars(dir, results_ids, star_descriptions: List[StarDescription]):
    newname = 'vsx_stars.txt'
    logging.info(f"Writing {newname}...")
    total_found = 0
    with open(dir+newname, 'wt') as fp:
        for vsx_id in results_ids:
            found = False if star_descriptions[vsx_id].path is '' else True
            total_found += 1 if found else 0
            fp.write(f"{vsx_id}{'' if found else '*'}:\t{star_descriptions[vsx_id].aavso_id}\t{utils.get_hms_dms(star_descriptions[vsx_id].coords)}\n")
        fp.write(f"# Total entries: {len(results_ids)}, found: {total_found}, not found: {len(results_ids) - total_found}\n")

def construct_star_descriptions(vastdir, wcs, all_stardict, selected_files, args):
    # start with the list of all detected stars
    all_stars = all_stardict.keys()

    # Start with the list of all measured stars
    measured_stars_dict = {}
    selected_files.sort()
    for afile in selected_files:
        star_id = get_starid_from_outfile(afile)
        measured_stars_dict[star_id] = (star_id, afile)
    #logging.debug(f"Star id list is : {measured_stars}")

    star_descriptions = do_calibration.get_empty_star_descriptions(all_stars)
    for idx, sd in enumerate(star_descriptions):
        #print(f"idx:{idx}, sd.localid: {sd.local_id}")
        sd.path = '' if sd.local_id not in measured_stars_dict else measured_stars_dict[sd.local_id][1]
        sd.xpos = all_stardict[int(sd.local_id)][0]
        sd.ypos = all_stardict[int(sd.local_id)][1]
        world_coords = wcs.all_pix2world(float(sd.xpos), float(sd.ypos), 0, ra_dec_order=True)
        logging.debug(f"world coords for star {sd.local_id}, {world_coords}")
        sd.coords = SkyCoord(world_coords[0], world_coords[1], unit='deg')

    if args.upsilon:
        # now we generate a list of StarDescriptions, with which all further processing will be done
        logging.info("Setting star_descriptions to upsilon candidates")
        star_descriptions = do_calibration.add_candidates_to_star_descriptions(star_descriptions, 0.1)

    if args.vsx:
        star_descriptions, results_ids = do_calibration.add_vsx_names_to_star_descriptions(star_descriptions, vsxcatalogdir, 0.01)
        logging.info(f"Added {len(results_ids)} vsx names to star descriptions")
        do_calibration.add_selected_match_to_stars(star_descriptions, results_ids, one_based=False) # select star ids
        # write the vsx stars used into a file
        results_ids.sort()
        write_vsx_stars(vastdir, results_ids, star_descriptions)

    if args.starfile:
        with open(settings.basedir + args.starfile, 'r') as fp:
            lines = fp.readlines()
            starlist = [x.rstrip() for x in lines]
            logging.debug(f"The list of stars read from the starfile is: {starlist} ")
            starlist = [int(x) for x in filter(str.isdigit, starlist)]
            logging.debug(f"The list of stars read from the starfile is: {starlist} ")
            logging.info(f"Selecting {starlist} stars added by {args.starfile}")
            do_calibration.add_selected_match_to_stars(star_descriptions, starlist, one_based=False) # select star ids

    return star_descriptions


# filter a list of star descriptions on the presence of a catalog
def catalog_filter(star: StarDescription, catalog_name):
    if star.get_catalog(catalog_name) is not None:
        return True
    return False


def interact():
    import code
    code.InteractiveConsole(locals=dict(globals(), **locals())).interact()


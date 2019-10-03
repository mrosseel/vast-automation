from multiprocessing import cpu_count
import logging
import re
import os.path
import subprocess
from functools import partial
from collections import namedtuple
import do_calibration
import do_charts_vast
import do_aavso_report
import do_charts_field
import do_compstars
import utils
from utils import get_star_description_cache
from reading import trash_and_recreate_dir
from reading import file_selector
from star_description import StarDescription
from photometry_blob import PhotometryBlob
from init_loader import init, settings
from astropy.coordinates import SkyCoord
from typing import List, Dict, Tuple

vsx_catalog_name = "vsx_catalog.bin"
vsxcatalogdir = '/home/jovyan/work/' + vsx_catalog_name

def run_do_rest(args):

    vastdir = utils.add_trailing_slash(args.datadir)
    fieldchartsdir = vastdir + 'fieldcharts/'
    aavsodir = vastdir + 'aavso/'
    logging.info(f"Directory where all files will be read & written: {vastdir}")
    wcs_file = vastdir+'new-image.fits'
    reference_frame = extract_reference_frame(vastdir)
    first_frame = extract_first_frame(vastdir)
    logging.info(f"The reference frame is {reference_frame}")
    logging.info(f"The first frame is {first_frame}")
    logging.info(f"reference header is {wcs_file}")
    while not os.path.isfile(wcs_file):
        logging.info(f"Please provide the reference header '{wcs_file}', which is an astrometry.net plate-solve of {reference_frame} and press Enter to continue...")
        subprocess.call("read -t 10", shell=True, executable='/bin/bash')

    # get wcs model from the reference header. Used in writing world positions and field charts
    wcs = do_calibration.get_wcs(wcs_file)

    all_stardict = read_stardict(vastdir)
    selected_files = file_selector(the_dir=vastdir, match_pattern="*.dat")
    logging.info(f"Number of found lightcurves: {len(selected_files)}, number of identified stars: {len(all_stardict.keys())}")
    star_descriptions = construct_star_descriptions(vastdir, wcs, all_stardict, selected_files, args)
    stardict = get_star_description_cache(star_descriptions)
    print("star descriptoios", star_descriptions[:100])
    write_augmented_autocandidates(vastdir, stardict)
    write_augmented_all_stars(vastdir, stardict)
    do_charts = args.lightcurve
    if args.field:
        do_charts_field.run_standard_field_charts(star_descriptions, wcs, fieldchartsdir, wcs_file)
    selected_stars = []
    if args.allstars:
        do_charts_vast.run(star_descriptions, vastdir, 'phase_all/', 'chart_all/', do_charts=do_charts,
                           nr_threads=cpu_count())
    else:
        if args.candidates:
            candidate_ids = get_autocandidates(vastdir)
            logging.info(f"Selecting {len(candidate_ids)} candidates to plot")
            selected_stars = do_calibration.select_star_descriptions(candidate_ids, star_descriptions)
            do_charts_vast.run(selected_stars, vastdir, 'phase_candidates/', 'chart_candidates/', do_charts=do_charts,
                               nr_threads=cpu_count())
        if args.vsx:
            vsx_stars = get_vsx_stars(star_descriptions)
            selected_stars = selected_stars + vsx_stars
            logging.info(f"Selecting {len(vsx_stars)} vsx stars to plot")
            do_charts_vast.run(vsx_stars, vastdir, 'phase_vsx/', 'chart_vsx/', do_charts=do_charts,
                               nr_threads=cpu_count())

    if args.aavso:
        # star_descriptions_ucac4 = do_calibration.add_ucac4_to_star_descriptions(star_descriptions)
        logging.info(f"AAVSO Reporting with: {len(selected_stars)} stars")
        print("avaso args", args.aavso)
        comparison_stars_1, comparison_stars_1_desc = do_compstars.get_fixed_compstars(star_descriptions, args.aavso)
        trash_and_recreate_dir(aavsodir)
        print("selected stars", selected_stars)
        # TODO put this in settings.txt
        sitelat = '-22 57 10'
        sitelong = '-68 10 49'
        sitealt= 2500
        for star in selected_stars:
            do_aavso_report.report(aavsodir, vastdir, star, sitelat, sitelong, sitealt,
                                   comparison_stars_1_desc[0], filter='V') # put filter to None for autodetect


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


# quickly test a few xy2sky conversions using our wcs and astropy
def wcs_test_pattern(wcs):
    logging.info("Outputting wcs testing pattern")
    test_pattern = [(0, 0), (500, 500), (136, 985), (-50, 500)]
    for tuple in test_pattern:
        logging.info(f"test pattern: {tuple}")
        result = wcs.all_pix2world(tuple[0], tuple[1], 0, ra_dec_order=True)
        logging.info(f"result: {result[0]}, {result[1]}")


def extract_reference_frame(from_dir):
    return extract_frame_helper(from_dir, "Ref.  image")


def extract_first_frame(from_dir):
    return extract_frame_helper(from_dir, "First image")


def extract_frame_helper(from_dir, marker):
    # Ref.  image: 2458586.50154 13.04.2019 00:00:41   ../../inputfiles/TXCar/fits/TXCar#45V_000601040_FLAT.fit
    result = [re.findall(marker + r': (?:.*) (.*)',line) for line in open(from_dir+'vast_summary.log')]
    return [x for x in result if x != []][0][0]


# get all possible stars with their x/y position from a log file
def read_stardict(vastdir):
    stardict = {}
    PixelPos = namedtuple('PixelPos', 'x y afile')
    for line in open(vastdir+'vast_list_of_all_stars.log'):
        splitline = line.split()
        stardict[int(splitline[0])] = PixelPos(splitline[1], splitline[2], f"out{splitline[0]}.dat")
    return stardict


# Note: this file seems to give incorrect xy positions wrt reference frame
# get all possible stars with their x/y position from a log file
# 14.460155 0.031190   215.230    19.626 out00007.dat
def read_data_m_sigma(vastdir) -> Dict[int, Tuple[int, int]]:
    stardict = {}
    PixelPos = namedtuple('PixelPos', 'x y afile')
    for line in open(vastdir+'data.m_sigma'):
        splitline = line.split()
        star_id = get_starid_from_outfile(splitline[4])
        stardict[star_id] = PixelPos(float(splitline[2]), float(splitline[3]), splitline[4])
    return stardict


def get_starid_from_outfile(outfile):
    m = re.search('out(.*).dat', outfile)
    return int(m.group(1).lstrip('0'))

    origname = 'vast_autocandidates.log'
    newname = 'vast_autocandidates_pos.txt'
    logging.info(f"Writing {newname}...")
    with open(dir+origname, 'r', encoding='utf-8') as infile, open(dir+newname, 'w') as outfile:
        for line in infile:
            linetext = line.rstrip()
            star_id = get_starid_from_outfile(linetext)
            if star_id in stardict:
                cacheentry = stardict[star_id]
                outfile.write(f"{linetext}{'' if cacheentry.path is not '' else '*'}\t{cacheentry.aavso_id}\t{utils.get_hms_dms(cacheentry.coords)}\n")
            else:
                outfile.write(f"{linetext}*\t{'None'}\n")


def get_autocandidates(dir) -> List[int]:
    origname = 'vast_autocandidates.log'
    result = []
    with open(dir+origname, 'r', encoding='utf-8') as infile:
        for line in infile:
            linetext = line.rstrip()
            star_id = get_starid_from_outfile(linetext)
            result.append(star_id)
    return result

def get_vsx_stars(stars: List[StarDescription]) -> List[StarDescription]:
    result = []
    for star in stars:
        if star.has_catalog("VSX"):
            result.append(star)
    return result

def write_augmented_autocandidates(dir, stardict: Dict[int, StarDescription]):
    origname = 'vast_autocandidates.log'
    newname = 'vast_autocandidates_pos.txt'
    logging.info(f"Writing {newname}...")
    with open(dir+origname, 'r', encoding='utf-8') as infile, open(dir+newname, 'w') as outfile:
        for line in infile:
            linetext = line.rstrip()
            star_id = get_starid_from_outfile(linetext)
            if star_id in stardict:
                cacheentry = stardict[star_id]
                outfile.write(f"{linetext}{'' if cacheentry.path is not '' else '*'}\t{cacheentry.aavso_id}\t{utils.get_hms_dms(cacheentry.coords)}\n")
            else:
                outfile.write(f"{linetext}*\t{'None'}\n")

def write_augmented_all_stars(dir, stardict: Dict[int, StarDescription]):
    origname = 'vast_list_of_all_stars.log'
    newname = 'vast_list_of_all_stars_pos.txt'
    logging.info(f"Writing {newname}...")
    with open(dir+origname, 'r', encoding='utf-8') as infile, open(dir+newname, 'w') as outfile:
        for line in infile:
            star_id = line.split()[0]
            if int(star_id) in stardict:
                cacheentry = stardict[int(star_id)]
                outfile.write(f"{star_id}\t{cacheentry.aavso_id}\t{utils.get_hms_dms(cacheentry.coords)}\t{cacheentry.coords.ra} {cacheentry.coords.dec}\n")

def write_vsx_stars(dir, results_ids, stars: List[StarDescription]):
    newname = 'vsx_stars.txt'
    logging.info(f"Writing {newname}...")
    total_found = 0
    stardict = utils.get_star_description_cache(stars)
    logging.debug(f"Receiving {len(stardict.keys())} as vsx input")
    with open(dir+newname, 'wt') as fp:
        for vsx_id in results_ids:
            current_sd = stardict[vsx_id]
            found = False if current_sd.path is '' else True
            assert vsx_id == current_sd.local_id
            total_found += 1 if found else 0
            fp.write(f"{vsx_id}{'' if found else '*'}:\t{current_sd.aavso_id}\t{utils.get_hms_dms(current_sd.coords)}\n")
        fp.write(f"# Total entries: {len(results_ids)}, found: {total_found}, not found: {len(results_ids) - total_found}\n")

def count_dat_entries(afile):
    return sum(1 for line in open(afile, 'r') if line.rstrip())

def construct_star_descriptions(vastdir, wcs, all_stardict, selected_files, args):

    # Start with the list of all measured stars
    measured_stars_dict = {}
    selected_files.sort()
    for afile in selected_files:
        star_id = get_starid_from_outfile(afile)
        measured_stars_dict[star_id] = afile
    #logging.debug(f"Star id list is : {measured_stars}")

    # intersect dicts
    intersect_dict = {x:(*all_stardict[x], measured_stars_dict[x]) for x in all_stardict if x in measured_stars_dict}
    logging.info(f"Calculating the intersect between all stars and measured stars, result has {len(intersect_dict)} entries.")
    # print("measuerd stars", measured_stars_dict)
    # print("all star dict", all_stardict.keys())

    star_descriptions = do_calibration.get_empty_star_descriptions(intersect_dict)
    for sd in star_descriptions:
        sd.path = '' if sd.local_id not in intersect_dict else intersect_dict[sd.local_id][3]
        sd.xpos = intersect_dict[int(sd.local_id)][0]
        sd.ypos = intersect_dict[int(sd.local_id)][1]
        sd.obs = -1
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
        print("vsx ids:", results_ids)
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

def has_option(obj, attr_name):
    return hasattr(obj, attr_name) and obj[attr_name]

# filter a list of star descriptions on the presence of a catalog
def catalog_filter(star: StarDescription, catalog_name):
    if star.get_catalog(catalog_name) is not None:
        return True
    return False


def interact():
    import code
    code.InteractiveConsole(locals=dict(globals(), **locals())).interact()


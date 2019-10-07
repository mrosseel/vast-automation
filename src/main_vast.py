from multiprocessing import cpu_count
import multiprocessing as mp
import tqdm
import logging
import re
import os.path
import subprocess
import numpy as np
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
from astropy.wcs import WCS
from typing import List, Dict, Tuple
from comparison_stars import ComparisonStars
from star_description import CatalogMatch
vsx_catalog_name = "vsx_catalog.bin"
vsxcatalogdir = '/home/jovyan/work/' + vsx_catalog_name
# star id -> xpos, ypos, filename
StarPosDict = Dict[str, Tuple[float, float, str]]

def run_do_rest(args):
    thread_count = cpu_count()-1
    vastdir = utils.add_trailing_slash(args.datadir)
    resultdir = utils.add_trailing_slash(args.resultdir) if args.resultdir is not None else vastdir
    fieldchartsdir = resultdir + 'fieldcharts/'
    aavsodir = resultdir + 'aavso/'
    logging.info(f"Directory where all files will be read: '{vastdir}' and written: '{resultdir}'")
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
    list_of_dat_files = file_selector(the_dir=vastdir, match_pattern="*.dat")
    logging.info(f"Number of found lightcurves: {len(list_of_dat_files)}, number of identified stars: {len(all_stardict.keys())}")
    star_descriptions = construct_star_descriptions(vastdir, resultdir, wcs, all_stardict, list_of_dat_files, args)
    stardict = get_star_description_cache(star_descriptions)
    logging.debug("First 10 star descriptions", star_descriptions[:10])
    write_augmented_autocandidates(vastdir, resultdir, stardict)
    write_augmented_all_stars(vastdir, resultdir, stardict)
    do_charts = args.lightcurve
    if args.field:
        do_charts_field.run_standard_field_charts(star_descriptions, wcs, fieldchartsdir, wcs_file)

    # load comparison stars
    comparison_stars_1, comparison_stars_1_desc = do_compstars.get_fixed_compstars(star_descriptions, args.aavso)
    comp_observations = []
    logging.info(f"Using comparison star ids:{comparison_stars_1}")
    for star in comparison_stars_1:
        read_comp_magdict = read_magdict(vastdir, star)
        # logging.info(f"Read comp magdict for {star}: {read_comp_magdict}")
        comp_observations.append(read_comp_magdict)
    comp_catalogmags = []
    comp_catalogerr = []
    for star in comparison_stars_1_desc:
        comp_catalogmags.append(star.vmag)
        comp_catalogerr.append(star.e_vmag)

    comp_stars = ComparisonStars(comparison_stars_1, comparison_stars_1_desc, comp_observations, comp_catalogmags, comp_catalogerr)
    logging.info(f"Comparison stars have {np.shape(comp_observations)}, {len(comp_observations[0])}, {len(comp_observations[1])} observations")
    candidate_stars = []
    if args.allstars:
        do_charts_vast.run(star_descriptions, comp_stars, vastdir, resultdir+'phase_all/', resultdir+'chart_all/', do_charts=do_charts,
                           nr_threads=thread_count)
        candidate_stars = star_descriptions
    else:
        if args.candidates:
            candidate_stars = do_calibration.get_catalog_stars(star_descriptions, "CANDIDATE")
            logging.info(f"Plotting {len(candidate_stars)} candidates...")
            do_charts_vast.run(candidate_stars, comp_stars, vastdir, resultdir+'phase_candidates/',
                               resultdir+'chart_candidates/', do_charts=do_charts, nr_threads=thread_count)
        if args.vsx:
            vsx_stars = do_calibration.get_catalog_stars(star_descriptions, "VSX")
            do_calibration.add_catalog_to_star_descriptions(vsx_stars, "SELECTED")
            logging.info(f"Plotting {len(vsx_stars)} vsx stars...")
            do_charts_vast.run(vsx_stars, comp_stars, vastdir, resultdir+'phase_vsx/', resultdir+'chart_vsx/', do_charts=do_charts,
                               nr_threads=thread_count)
        if args.starfile:
            selected_stars = do_calibration.get_catalog_stars(star_descriptions, "STARFILE")
            print("Selected stars from starfile", selected_stars)
            do_charts_vast.run(selected_stars, comp_stars, vastdir, resultdir+'phase_selected/',
                               resultdir+'chart_selected/', do_charts=do_charts, nr_threads=thread_count)

    if args.aavso:
        # star_descriptions_ucac4 = do_calibration.add_ucac4_to_star_descriptions(star_descriptions)
        logging.info(f"AAVSO Reporting with: {len(selected_stars)} stars and comparison stars {args.aavso}")

        trash_and_recreate_dir(aavsodir)
        # TODO put this in settings.txt
        sitelat = '-22 57 10'
        sitelong = '-68 10 49'
        sitealt = 2500
        observer = 'RMH'
        pool = mp.Pool(thread_count, maxtasksperchild=10)
        # put filter to None for autodetect
        func = partial(do_aavso_report.report, target_dir=aavsodir, vastdir=vastdir,
                       sitelat=sitelat, sitelong=sitelong, sitealt=sitealt,
                       comparison_star=comparison_stars_1_desc[0], filter='V', observer=observer, chunk_size=args.aavsolimit)
        for _ in tqdm.tqdm(pool.imap_unordered(func, selected_stars, 5), total=len(selected_stars), unit="files"):
            pass


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


# get a dict with star_id -> xpos ypos filename
def read_stardict(vastdir: str) -> StarPosDict:
    stardict = {}
    PixelPos = namedtuple('PixelPos', 'x y afile')
    for line in open(vastdir+'vast_list_of_all_stars.log'):
        splitline = line.split()
        stardict[int(splitline[0])] = PixelPos(splitline[1], splitline[2], f"out{splitline[0]}.dat")
    return stardict


# get all possible stars with dict: JD, (mag, error)
def read_magdict(vastdir, star_id):
    stardict = {}
    starfile = f"{vastdir}{star_to_dat(star_id)}"
    for line in open(starfile):
        splitline = line.split()
        # {JD, (mag, magerr)}
        stardict[str(splitline[0])] = (splitline[1], splitline[2])
    return stardict


def star_to_dat(star: int):
    return f"out{star:05}.dat"

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


def get_autocandidates(dir) -> List[int]:
    origname = 'vast_autocandidates.log'
    result = []
    with open(dir+origname, 'r', encoding='utf-8') as infile:
        for line in infile:
            linetext = line.rstrip()
            star_id = get_starid_from_outfile(linetext)
            result.append(star_id)
    return result


def write_augmented_autocandidates(readdir: str, writedir: str, stardict: Dict[int, StarDescription]):
    origname = f"{readdir}vast_autocandidates.log"
    newname = f"{writedir}vast_autocandidates_pos.txt"
    logging.info(f"Writing {newname}...")
    with open(origname, 'r', encoding='utf-8') as infile, open(newname, 'w') as outfile:
        for line in infile:
            linetext = line.rstrip()
            star_id = get_starid_from_outfile(linetext)
            if star_id in stardict:
                cacheentry = stardict[star_id]
                outfile.write(f"{linetext}{'' if cacheentry.path is not '' else '*'}\t{cacheentry.aavso_id}\t{utils.get_lesve_coords(cacheentry.coords)}\n")
            else:
                outfile.write(f"{linetext}*\t{'None'}\n")


def write_augmented_all_stars(readdir: str, writedir: str, stardict: Dict[int, StarDescription]):
    origname = f"{readdir}vast_list_of_all_stars.log"
    newname = f"{writedir}vast_list_of_all_stars_pos.txt"
    logging.info(f"Writing {newname}...")
    with open(origname, 'r', encoding='utf-8') as infile, open(newname, 'w') as outfile:
        for line in infile:
            star_id = line.split()[0]
            if int(star_id) in stardict:
                cacheentry = stardict[int(star_id)]
                outfile.write(f"{star_id}\t{cacheentry.aavso_id}\t{utils.get_hms_dms(cacheentry.coords)}\t{cacheentry.coords.ra} {cacheentry.coords.dec}\n")


def write_vsx_stars(resultdir, results_ids, stars: List[StarDescription]):
    newname = f"{resultdir}vsx_stars.txt"
    logging.info(f"Writing {newname}...")
    total_found = 0
    stardict = utils.get_star_description_cache(stars)
    logging.debug(f"Receiving {len(stardict.keys())} as vsx input")
    with open(newname, 'wt') as fp:
        for vsx_id in results_ids:
            current_sd = stardict[vsx_id]
            found = False if current_sd.path is '' else True
            assert vsx_id == current_sd.local_id
            total_found += 1 if found else 0
            fp.write(f"{vsx_id}{'' if found else '*'}:\t{current_sd.aavso_id}\t{utils.get_lesve_coords(current_sd.coords)}\n")
        fp.write(f"# Total entries: {len(results_ids)}, found: {total_found}, not found: {len(results_ids) - total_found}\n")


def count_dat_entries(afile):
    return sum(1 for line in open(afile, 'r') if line.rstrip())


# constructs a list of star descriptions with catalog matches according to args
def construct_star_descriptions(vastdir: str, resultdir: str, wcs: WCS, all_stardict: StarPosDict, list_of_dat_files: List[str], args):
    # Start with the list of all measured stars
    stars_with_file_dict = {}
    list_of_dat_files.sort()
    for afile in list_of_dat_files:
        star_id = get_starid_from_outfile(afile)
        stars_with_file_dict[star_id] = afile

    # intersect dict, results in starid -> (xpos, ypos, shortfile, longfile),
    # example:  42445: ('175.948', '1194.074', 'out42445.dat', 'support/vast-1.0rc84/out42445.dat')
    intersect_dict = {x:(*all_stardict[x], stars_with_file_dict[x]) for x in all_stardict if x in stars_with_file_dict}
    logging.info(f"Calculating the intersect between all stars and measured stars, result has {len(intersect_dict)} entries.")

    # get SD's for all stars which are backed by a file with measurements
    star_descriptions = do_calibration.get_empty_star_descriptions(intersect_dict)
    for sd in star_descriptions:
        sd.path = '' if sd.local_id not in intersect_dict else intersect_dict[sd.local_id][3]
        sd.xpos = intersect_dict[int(sd.local_id)][0]
        sd.ypos = intersect_dict[int(sd.local_id)][1]
        sd.obs = -1
        world_coords = wcs.all_pix2world(float(sd.xpos), float(sd.ypos), 0, ra_dec_order=True)
        logging.debug(f"world coords for star {sd.local_id}, {world_coords}")
        sd.coords = SkyCoord(world_coords[0], world_coords[1], unit='deg')

    # Add VSX information to SDs
    star_descriptions, results_ids = do_calibration.add_vsx_names_to_star_descriptions(star_descriptions,
                                                                                       vsxcatalogdir, 0.01)
    logging.info(f"Added {len(results_ids)} vsx names to star descriptions")

    # write the vsx stars used into a file
    results_ids.sort()
    write_vsx_stars(resultdir, results_ids, star_descriptions)

    if args.candidates:
        candidate_ids = get_autocandidates(vastdir)
        logging.info(f"Tagging {len(candidate_ids)} candidates")
        selected_stars = do_calibration.select_star_descriptions(candidate_ids, star_descriptions)
        for star in selected_stars:
            star.match.append(CatalogMatch(name_of_catalog="CANDIDATE", catalog_id=star.local_id,
                                           name=star.local_id))
        do_calibration.add_catalog_to_star_descriptions(selected_stars, "SELECTED")

    if args.starfile:
        with open(args.starfile, 'r') as fp:
            lines = fp.readlines()
            starlist = [x.rstrip() for x in lines]
            logging.info(f"The list of stars read from the starfile is: {starlist} ")
            starlist = [int(x) for x in filter(str.isdigit, starlist)]
            logging.info(f"The list of stars read from the starfile is: {starlist} ")
            logging.info(f"Selecting {starlist} stars added by {args.starfile}")
            starfile_stars = list(filter(lambda x: x.local_id in starlist, star_descriptions))
            print("starfile stars", starfile_stars)
            do_calibration.add_catalog_to_star_descriptions(starfile_stars, "SELECTED")
            do_calibration.add_catalog_to_star_descriptions(starfile_stars, "STARFILE")
            print("starfile stars", starfile_stars)
    return star_descriptions


def has_option(obj, attr_name):
    return hasattr(obj, attr_name) and obj[attr_name]


def interact():
    import code
    code.InteractiveConsole(locals=dict(globals(), **locals())).interact()


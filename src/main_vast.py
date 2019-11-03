from multiprocessing import cpu_count
import tqdm
import logging
import re
import os
import os.path
import subprocess
import numpy as np
from functools import partial
from collections import namedtuple
import do_calibration
import do_charts_vast
import do_charts_field
import do_compstars
import reading
import star_metadata
import utils
import star_description
from utils import get_star_description_cache
from reading import trash_and_recreate_dir
from reading import file_selector
from star_description import StarDescription
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.wcs import WCS
from typing import List, Dict, Tuple
from comparison_stars import ComparisonStars
from pathlib import PurePath
from ucac4 import UCAC4
import hugo_site
import pandas as pd
from star_metadata import SelectedStarData, CatalogData, StarFileData, CompStarData

vsx_catalog_name = "vsx_catalog.bin"
vsxcatalogdir = PurePath(os.getcwd(), vsx_catalog_name)
# star id -> xpos, ypos, filename
StarPosDict = Dict[str, Tuple[float, float, str]]
StarDict = Dict[int, StarDescription]
STAR_KEEPER_PERCENTAGE = 0.1


def run_do_rest(args):
    thread_count = cpu_count() - 1
    vastdir = utils.add_trailing_slash(args.datadir)
    resultdir = clean_and_create_resultdir(args.resultdir, vastdir)
    fieldchartsdir = resultdir + 'fieldcharts/'
    aavsodir = resultdir + 'aavso/'
    do_charts = args.light
    do_phase = args.phase
    do_aavso = args.aavso
    logging.info(f"Directory where all files will be read: '{vastdir}' and written: '{resultdir}'")
    wcs_file = vastdir + 'new-image.fits'
    reference_frame = extract_reference_frame(vastdir)
    first_frame = extract_first_frame(vastdir)
    frames_used = int(extract_images_used(vastdir))
    logging.info(f"{frames_used} frames were used for photometry")
    logging.info(f"The reference frame is {reference_frame}")
    logging.info(f"The first frame is {first_frame}")
    logging.info(f"reference header is {wcs_file}")
    while not os.path.isfile(wcs_file):
        logging.info(
            f"Please provide the reference header '{wcs_file}', which is an astrometry.net plate-solve of {reference_frame} and press Enter to continue...")
        subprocess.call("read -t 10", shell=True, executable='/bin/bash')

    # get wcs model from the reference header. Used in writing world positions and field charts
    wcs = do_calibration.get_wcs(wcs_file)

    all_stardict = read_stardict(vastdir)
    list_of_dat_files = file_selector(the_dir=vastdir, match_pattern="*.dat")
    logging.info(
        f"Number of found lightcurves: {len(list_of_dat_files)}, number of identified stars: {len(all_stardict.keys())}")
    star_descriptions = construct_star_descriptions(vastdir, resultdir, wcs, all_stardict, list_of_dat_files,
                                                    frames_used, args)
    stardict = get_star_description_cache(star_descriptions)
    logging.debug("First (max) 10 star descriptions",
                  star_descriptions[:10] if len(star_descriptions) >= 10 else star_descriptions)
    write_augmented_autocandidates(vastdir, resultdir, stardict)
    write_augmented_all_stars(vastdir, resultdir, stardict)
    candidate_stars = utils.get_stars_with_metadata(star_descriptions, "CANDIDATE", exclude="VSX")
    rmhhmb_stars = utils.get_stars_with_metadata(star_descriptions, "RMH+HMB")
    candidate_stars = utils.add_star_lists(candidate_stars, rmhhmb_stars)
    vsx_stars = utils.get_stars_with_metadata(star_descriptions, "VSX")
    starfile_stars = utils.get_stars_with_metadata(star_descriptions, "STARFILE")
    write_augmented_starfile(resultdir, starfile_stars)
    interesting_stars = starfile_stars + vsx_stars + candidate_stars

    # set ucac4 stars selected
    ucac4 = UCAC4()
    ucac4.add_ucac4_to_sd(starfile_stars)

    comp_stars = read_comparison_stars(star_descriptions, args.checkstarfile, vastdir, stardict, ucac4)

    # Set comp stars for every star
    logging.info("Setting per star comparison stars...")
    if args.checkstarfile:
        utils.add_metadata(star_descriptions, CompStarData(compstar_ids=comp_stars.ids))
    else:
        do_compstars.add_closest_compstars(interesting_stars, comp_stars, 10)

    if args.allstars:
        do_charts_vast.run(star_descriptions, comp_stars, vastdir, resultdir, 'phase_all/', 'light_all/', 'aavso_all/',
                           do_phase=do_phase, do_charts=do_charts, do_aavso=do_aavso, nr_threads=thread_count,
                           desc="Phase/light/aavso of ALL stars")
    else:
        if args.candidates:
            logging.info(f"Plotting {len(candidate_stars)} candidates...")
            do_charts_vast.run(candidate_stars, comp_stars, vastdir, resultdir, 'phase_candidates/',
                               'light_candidates/', 'aavso_candidates/', do_phase=do_phase, do_charts=do_charts,
                               do_aavso=do_aavso, nr_threads=thread_count, desc="Phase/light/aavso of candidates")
        if args.vsx:
            do_calibration.add_metadata_to_star_descriptions(vsx_stars, ["SELECTED"], strict=False)
            logging.info(f"Plotting {len(vsx_stars)} vsx stars...")
            do_charts_vast.run(vsx_stars, comp_stars, vastdir, resultdir, 'phase_vsx/', 'light_vsx/', 'aavso_vsx/',
                               do_phase=do_phase, do_charts=do_charts, do_aavso=do_aavso, nr_threads=thread_count,
                               desc="Phase/light/aavso of VSX stars")
        if args.selectedstarfile:
            do_charts_vast.run(starfile_stars, comp_stars, vastdir, resultdir, 'phase_selected/', 'light_selected/',
                               'aavso_selected', do_phase=do_phase, do_charts=do_charts,
                               do_aavso=do_aavso, nr_threads=thread_count, desc="Phase/light/aavso of selected stars")

    if args.field:
        do_charts_field.run_standard_field_charts(star_descriptions, wcs, fieldchartsdir, wcs_file, comp_stars)

    if args.site:
        ids = [x.local_id for x in starfile_stars]
        logging.info(f"Creating site entry with these {len(starfile_stars)} selected stars: {ids}")
        hugo_site.run(args.site, starfile_stars, resultdir)


def read_comparison_stars(star_descriptions: List[StarDescription], checkstarfile: str, vastdir: str,
                          stardict: StarDict, ucac4: UCAC4) -> ComparisonStars:
    if checkstarfile:
        # load comparison stars
        checkstars = read_checkstars(checkstarfile)
        comparison_stars_1, comparison_stars_1_desc = do_compstars.get_fixed_compstars(star_descriptions, checkstars)
    else:
        ucac4.add_ucac4_to_sd(star_descriptions)
        comparison_stars_1, comparison_stars_1_desc = do_compstars.get_calculated_compstars(vastdir, stardict)
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
    comp_stars = ComparisonStars(comparison_stars_1, comparison_stars_1_desc, comp_observations, comp_catalogmags,
                                 comp_catalogerr)
    logging.info(
        f"Comparison stars have {np.shape(comp_observations)}, {len(comp_observations[0])}, {len(comp_observations[1])} observations")
    return comp_stars


def clean_and_create_resultdir(argsdir: str, vastdir: str):
    resultdir = utils.add_trailing_slash(argsdir) if argsdir is not None else vastdir
    # if resultdir does not exist, create it
    if not os.path.isdir(resultdir):
        logging.info(f"The resultdir '{resultdir}' does not exist, creating it...")
        reading.create_dir(resultdir)
    return resultdir


# quickly test a few xy2sky conversions using our wcs and astropy
def wcs_test_pattern(wcs):
    logging.info("Outputting wcs testing pattern")
    test_pattern = [(0, 0), (500, 500), (136, 985), (-50, 500)]
    for tuple in test_pattern:
        logging.info(f"test pattern: {tuple}")
        result = wcs.all_pix2world(tuple[0], tuple[1], 0, ra_dec_order=True)
        logging.info(f"result: {result[0]}, {result[1]}")


def extract_reference_frame(from_dir):
    return extract_frame_from_summary_helper(from_dir, "Ref.  image")


def extract_first_frame(from_dir):
    return extract_frame_from_summary_helper(from_dir, "First image")


def extract_images_used(from_dir):
    result = [re.findall(r'Images used for photometry (.*)', line) for line in open(from_dir + 'vast_summary.log')]
    return [x for x in result if x != []][0][0]


def extract_frame_from_summary_helper(from_dir, marker):
    # Ref.  image: 2458586.50154 13.04.2019 00:00:41   ../../inputfiles/TXCar/fits/TXCar#45V_000601040_FLAT.fit
    result = [re.findall(marker + r': (?:.*) (.*)', line) for line in open(from_dir + 'vast_summary.log')]
    return [x for x in result if x != []][0][0]


# get a dict with star_id -> xpos ypos filename
def read_stardict(vastdir: str) -> StarPosDict:
    stardict = {}
    PixelPos = namedtuple('PixelPos', 'x y afile')
    for line in open(vastdir + 'vast_list_of_all_stars.log'):
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
        stardict[str(splitline[0])] = (float(splitline[1]), float(splitline[2]))
    return stardict


def star_to_dat(star: int):
    return f"out{star:05}.dat"


# Note: this file seems to give incorrect xy positions wrt reference frame
# get all possible stars with their x/y position from a log file
# 14.460155 0.031190   215.230    19.626 out00007.dat
def read_data_m_sigma(vastdir) -> Dict[int, Tuple[int, int]]:
    stardict = {}
    PixelPos = namedtuple('PixelPos', 'x y afile')
    for line in open(vastdir + 'data.m_sigma'):
        splitline = line.split()
        star_id = utils.get_starid_from_outfile(splitline[4])
        stardict[star_id] = PixelPos(float(splitline[2]), float(splitline[3]), splitline[4])
    return stardict


def read_checkstars(checkstar_file: str) -> List[str]:
    result = []
    for line in open(checkstar_file):
        result.append(line.strip())
    return result


def get_autocandidates(dir: str) -> List[int]:
    origname = 'vast_autocandidates.log'
    result = []
    with open(PurePath(dir, origname), 'r', encoding='utf-8') as infile:
        for line in infile:
            linetext = line.rstrip()
            star_id = utils.get_starid_from_outfile(linetext)
            result.append(star_id)
    return result


def write_augmented_autocandidates(readdir: str, writedir: str, stardict: StarDict):
    origname = f"{readdir}vast_autocandidates.log"
    newname = f"{writedir}vast_autocandidates_pos.txt"
    logging.info(f"Writing {newname}...")
    with open(origname, 'r', encoding='utf-8') as infile, open(newname, 'w') as outfile:
        for line in infile:
            linetext = line.rstrip()
            star_id = utils.get_starid_from_outfile(linetext)
            if star_id in stardict:
                cacheentry = stardict[star_id]
                outfile.write(
                    f"{linetext}{'' if cacheentry.path is not '' else '*'}\t{cacheentry.aavso_id}\t{utils.get_lesve_coords(cacheentry.coords)}\n")
            else:
                outfile.write(f"{linetext}*\t{'None'}\n")


def write_augmented_all_stars(readdir: str, writedir: str, stardict: StarDict):
    origname = f"{readdir}vast_list_of_all_stars.log"
    newname = f"{writedir}vast_list_of_all_stars_pos.txt"
    logging.info(f"Writing {newname}...")
    with open(origname, 'r', encoding='utf-8') as infile, open(newname, 'w') as outfile:
        for line in infile:
            star_id = line.split()[0]
            if int(star_id) in stardict:
                cacheentry = stardict[int(star_id)]
                outfile.write(
                    f"{star_id}\t{cacheentry.aavso_id}\t{utils.get_hms_dms(cacheentry.coords)}\t{cacheentry.coords.ra} {cacheentry.coords.dec}\n")


# naam, ra, dec, max, min, type, periode, epoch?
def write_augmented_starfile(writedir: str, starfile_stars: List[StarDescription]):
    newname = f"{writedir}starfile.txt"
    logging.info(f"Writing {newname}...")
    sorted_stars = utils.sort_rmh_hmb(starfile_stars)
    with open(newname, 'w') as outfile:
        outfile.write(f"# our_name,ra,dec,minmax,var_type,period,period_err,epoch\n")
        for star in sorted_stars:
            metadata: StarFileData = star.get_metadata("STARFILE")
            outfile.write(
                f"{metadata.our_name},{star.coords.ra.deg:.5f},{star.coords.dec.deg:.5f},{metadata.minmax},"
                f"{metadata.var_type},{metadata.period},{metadata.period_err},{metadata.epoch}\n")


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
            fp.write(
                f"{vsx_id}{'' if found else '*'}:\t{current_sd.aavso_id}\t{utils.get_lesve_coords(current_sd.coords)}\n")
        fp.write(
            f"# Total entries: {len(results_ids)}, found: {total_found}, not found: {len(results_ids) - total_found}\n")


def count_dat_entries(afile):
    return sum(1 for line in open(afile, 'r') if line.rstrip())


# constructs a list of star descriptions with catalog matches according to args
def construct_star_descriptions(vastdir: str, resultdir: str, wcs: WCS, all_stardict: StarPosDict,
                                list_of_dat_files: List[str], frames_used: int, args):
    # Start with the list of all measured stars
    stars_with_file_dict = {}
    list_of_dat_files.sort()
    for afile in list_of_dat_files:
        star_id = utils.get_starid_from_outfile(afile)
        stars_with_file_dict[star_id] = afile

    # intersect dict, results in starid -> (xpos, ypos, shortfile, longfile),
    # example:  42445: ('175.948', '1194.074', 'out42445.dat', 'support/vast-1.0rc84/out42445.dat')
    intersect_dict = {x: (*all_stardict[x], stars_with_file_dict[x]) for x in all_stardict if x in stars_with_file_dict}
    logging.info(
        f"Calculating the intersect between all stars and measured stars, result has {len(intersect_dict)} entries.")

    # get SD's for all stars which are backed by a file with measurements
    star_descriptions = do_calibration.get_empty_star_descriptions(intersect_dict)
    for sd in star_descriptions:
        sd.path = '' if sd.local_id not in intersect_dict else intersect_dict[sd.local_id][3]
        sd.xpos = intersect_dict[int(sd.local_id)][0]
        sd.ypos = intersect_dict[int(sd.local_id)][1]
        sd.obs = -1
        world_coords = wcs.all_pix2world(float(sd.xpos), float(sd.ypos), 0, ra_dec_order=True)
        # logging.debug(f"world coords for star {sd.local_id}, {world_coords}")
        sd.coords = SkyCoord(world_coords[0], world_coords[1], unit='deg')

    # add line counts
    pool = utils.get_pool(cpu_count() * 2,
                          maxtasksperchild=1000)  # lot's of small files, needs many threads to fill cpu
    stars = []
    for star in tqdm.tqdm(pool.imap_unordered(set_lines, star_descriptions, 5), total=len(star_descriptions),
                          desc="Counting obs per star", unit="files"):
        stars.append(star)
    star_descriptions = stars
    # only keep stars which are present on at least 10% of the images
    star_descriptions = list(filter(lambda x: x.obs > frames_used * STAR_KEEPER_PERCENTAGE, star_descriptions))
    logging.info(f"Filtered star descriptions to {len(star_descriptions)} stars")
    stardict = get_star_description_cache(star_descriptions)

    # Add VSX information to SDs
    star_descriptions, results_ids = do_calibration.add_vsx_names_to_star_descriptions(star_descriptions,
                                                                                       vsxcatalogdir, 0.01)
    logging.info(f"Added {len(results_ids)} vsx names to star descriptions")
    test = utils.get_stars_with_metadata(star_descriptions, "VSX")
    logging.info(f"Test Tagged {len(test)} stars as VSX.")

    # write the vsx stars used into a file
    results_ids.sort()
    write_vsx_stars(resultdir, results_ids, star_descriptions)

    # tag all candidates with a 'selected' and 'candidate' catalog
    tag_candidates(vastdir, star_descriptions)

    if args.selectedstarfile:
        tag_starfile(args.selectedstarfile, stardict)
        logging.info(f"Succesfully read {len(list(filter(lambda x: x.has_metadata('STARFILE'), star_descriptions)))} "
                     f"stars from file:"
                     f" {[x.local_id for x in list(filter(lambda x: x.has_metadata('STARFILE'), star_descriptions))]}")

    if args.owncatalog:
        tag_owncatalog(args.owncatalog, star_descriptions)

    # add ucac4 id's
    starfile_stars = utils.get_stars_with_metadata(star_descriptions, "STARFILE")
    return star_descriptions


def tag_candidates(vastdir: str, star_descriptions: List[StarDescription]):
    candidate_ids = get_autocandidates(vastdir)
    selected_stars = do_calibration.select_star_descriptions(candidate_ids, star_descriptions)
    do_calibration.add_metadata_to_star_descriptions(selected_stars, ["SELECTED", "CANDIDATE"], strict=False)


def tag_starfile(selectedstarfile: str, stardict: StarDict):
    try:
        df = pd.read_csv(selectedstarfile, delimiter=',', comment='#',
                         names=['local_id', 'minmax', 'var_type', 'our_name', 'period', 'period_err', 'epoch'],
                         dtype={'local_id': int, 'minmax': str, 'period': float, 'period_err': float},
                         skipinitialspace=True)
        df = df.replace({np.nan: None})
        logging.info(f"Selecting {len(df)} stars added by {selectedstarfile}, {df['local_id'].to_numpy()}")
        for idx, row in df.iterrows():
            the_star = stardict.get(row['local_id'])
            if the_star is None:
                logging.error(f"Could not find star {row['local_id']}, consider removing it from your txt file")
                continue
            the_star.metadata = StarFileData(row['local_id'], row['minmax'], row['var_type'], row['our_name'],
                                             row['period'], row['period_err'], row['epoch'])
            the_star.metadata = SelectedStarData()
            logging.info(f"starfile {the_star.local_id} metadata: {the_star.metadata}, "
                         f"{the_star.get_metadata('STARFILE')}")
            logging.info(f"starfile {the_star.get_metadata('STARFILE')}")
        logging.info(f"Tagged {len(df)} stars as selected by file.")
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        logging.error(message)
        logging.error(f"Could not read {selectedstarfile}, star {row['local_id']}")


def tag_owncatalog(owncatalog: str, stars: List[StarDescription]):
    # outfile.write(f"# our_name,ra,dec,minmax,var_type,period,epoch\n")
    df = pd.read_csv(owncatalog, delimiter=',', comment='#',
                     names=['our_name', 'ra', 'dec', 'minmax', 'var_type', 'period', 'epoch'],
                     dtype={'ra': float, 'dec': float, 'minmax': str},
                     skipinitialspace=True)
    df = df.replace({np.nan: None})
    skycoord: SkyCoord = do_calibration.create_generic_astropy_catalog(df['ra'], df['dec'])
    star_catalog = do_calibration.create_star_descriptions_catalog(stars)
    idx, d2d, d3d = match_coordinates_sky(skycoord, star_catalog, nthneighbor=1)
    for count, index in enumerate(idx):
        entry = df.iloc[count]
        stars[index].metadata = CatalogData(key="RMH+HMB", catalog_id=entry['our_name'],
                                            name=entry['our_name'],
                                            coords=SkyCoord(entry['ra'], entry['dec'], unit="deg"),
                                            separation=d2d[count].degree)
        if d2d[count].degree > 0.01:
            logging.warning(f"Separation between {df.iloc(count)['our_name']} "
                            f"and {stars[index].local_id} is {d2d[count]}")


def set_lines(star: StarDescription):
    star.obs = sum(1 for line in open(star.path) if line.rstrip())
    return star


def has_option(obj, attr_name):
    return hasattr(obj, attr_name) and obj[attr_name]


def interact():
    import code
    code.InteractiveConsole(locals=dict(globals(), **locals())).interact()

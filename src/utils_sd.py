from typing import List

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from pathlib import Path

import reading
import utils
import logging

# returns list of star descriptions
from star_description import StarDescription


def get_empty_star_descriptions(star_id_list=None):
    # returns {'name': [ra.deg, dec.deg ]}
    result = []
    for star in star_id_list:
        result.append(StarDescription(local_id=star))
    return result


def construct_star_descriptions(vastdir: str, stars: List[int] = None):
    vastdir = utils.add_trailing_slash(vastdir)
    wcs_file, wcs = reading.read_wcs_file(vastdir)
    all_files = utils.file_selector(the_dir=vastdir, match_pattern="*.dat")
    if stars is not None:
        all_files = [vastdir + reading.star_to_dat(int(x)) for x in stars]
    sds = construct_raw_star_descriptions(vastdir, wcs, all_files)
    return sds


def construct_raw_star_descriptions(vastdir: str, wcs: WCS, list_of_dat_files: List[str] = None,
                                    star_keeper_percentage=0.1):
    all_stardict = reading.starid_to_xy_file_dict(vastdir)
    # if no list of dat files is passed, we use all the dat files
    if list_of_dat_files is None:
        list_of_dat_files = utils.file_selector(the_dir=vastdir, match_pattern="*.dat")
    logging.info(
        f"Number of found lightcurves: {len(list_of_dat_files)}, "
        f"number of identified stars: {len(all_stardict.keys())}")

    # Start with the list of all measured stars
    stars_with_file_dict = {}
    list_of_dat_files.sort()
    for afile in list_of_dat_files:
        star_id = utils.get_starid_from_outfile(afile)
        stars_with_file_dict[star_id] = afile

    # intersect dict, results in starid -> (xpos, ypos, shortfile, longfile),
    # example:  42445: ('175.948', '1194.074', 'out42445.dat', 'support/vast-1.0rc84/out42445.dat')
    intersect_dict = {x: (*all_stardict[x], stars_with_file_dict[x]) for x in all_stardict if
                      x in stars_with_file_dict}
    logging.info(
        f"Calculating the intersect between all stars and measured stars, result has {len(intersect_dict)} entries.")

    # get SD's for all stars which are backed by a file with measurements
    star_descriptions = get_empty_star_descriptions(intersect_dict)
    obsdict = reading.count_number_of_observations(vastdir)
    for sd in star_descriptions:
        sd.path = '' if sd.local_id not in intersect_dict else intersect_dict[sd.local_id][3]
        sd.xpos = intersect_dict[int(sd.local_id)][0]
        sd.ypos = intersect_dict[int(sd.local_id)][1]
        path_filename = Path(sd.path).name
        sd.obs = obsdict[path_filename] if path_filename in obsdict else -1
        world_coords = wcs.all_pix2world(float(sd.xpos), float(sd.ypos), 0, ra_dec_order=True)
        # logging.debug(f"world coords for star {sd.local_id}, {world_coords}")
        sd.coords = SkyCoord(world_coords[0], world_coords[1], unit='deg')

    frames_used = int(reading.extract_images_used(vastdir))
    logging.info(f"{frames_used} frames were used for photometry")

    # only keep stars which are present on at least 10% of the images
    star_descriptions = list(filter(lambda x: x.obs > frames_used * star_keeper_percentage, star_descriptions))
    return star_descriptions

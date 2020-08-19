import os
import reading
import logging
from star_description import StarDescription
from star_metadata import CatalogData
from star_metadata import UpsilonData
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy import units as u
from astropy.coordinates import Angle

# astroquery imports the keyring backend, but why?
# from astroquery.vizier import Vizier
import vsx_pickle
import do_calibration
import numpy as np
import pandas as pd
from functools import partial
from multiprocessing import Pool
from typing import List
import utils
from collections import namedtuple

from star_metadata import StarMetaData


def get_wcs(wcs_file):
    hdulist = fits.open(wcs_file)
    # data = hdulist[0].data.astype(float)
    header = hdulist[0].header
    wcs = WCS(header)
    return wcs


############# star description utils #################


def select_star_descriptions(star_id_list: List[int], stars: List[StarDescription]):
    return [x for x in stars if x.local_id in star_id_list]


# for testinsg
def get_random_star_descriptions(nr=10):
    result = []
    for idx in range(nr):
        result.append(
            StarDescription(local_id=idx, coords=SkyCoord(idx, idx, unit="deg"))
        )
    return result


# given a list of SD's and a list of star id's, output the matching SD's in log.debug
def log_star_descriptions(star_descriptions: StarDescription, star_id_list):
    logging.debug(f"Logging {len(star_id_list)} star  descriptions:")
    for star_id in star_id_list:
        logging.debug(star_descriptions[star_id])


############# star description utils #################

# adds UpsilonMatch to existing StarDescription
def add_candidates_to_star_descriptions(
    stars: List[StarDescription], threshold_prob=0.5, check_flag=False
):
    df = get_upsilon_candidates_raw(threshold_prob, check_flag)
    if df is None:
        return stars
    for index, row in df.iterrows():
        upsilon_match = UpsilonData(
            var_type=row["label"],
            probability=row["probability"],
            flag=row["flag"],
            period=row["period"],
        )
        stars[index].metadata.append(upsilon_match)
    return stars


# common upsilon code
def get_upsilon_candidates_raw(threshold_prob, check_flag):
    upsilon_file = settings.basedir + "upsilon_output.txt"
    try:
        df = pd.read_csv(upsilon_file, index_col=0)
    except:
        logging.error(
            f"could not read {upsilon_file}, skipping upsilon candidate processing"
        )
        return None
    df.sort_values(by="probability", ascending=False)
    df = df[df["label"] != "NonVar"]
    df = df[df["probability"] > threshold_prob]
    if check_flag:
        df = df[df["flag"] != 1]
    return df


# returns {'star_id': [label, probability, flag, period, SkyCoord, match_name, match_skycoord,
# match_type, separation_deg]}
def add_vsx_names_to_star_descriptions(
    star_descriptions: List[StarDescription], vsxcatalogdir: str, max_separation=0.01
):
    result_ids = []
    # copy.deepcopy(star_descriptions)
    vsx_catalog, vsx_dict = create_vsx_astropy_catalog(vsxcatalogdir)
    star_catalog = create_star_descriptions_catalog(star_descriptions)
    # vsx catalog is bigger in this case than star_catalog, but if we switch then different stars can be
    # matched with the same variable which is wrong.
    #
    # idx : integer array
    # Indices into catalogcoord to get the matched points for each matchcoord. Shape matches matchcoord.
    #
    # sep2d : Angle
    # The on-sky separation between the closest match for each matchcoord and the matchcoord. Shape matches matchcoord.
    idx, d2d, _ = match_coordinates_sky(star_catalog, vsx_catalog)
    logging.debug(f"length of idx: {len(idx)}")
    results_dict = {}  # have temp results_dict so we can remove duplicates
    VsxInfo = namedtuple("VsxInfo", "starid_0 vsx_id sep")
    for index_star_catalog, entry in enumerate(d2d):
        if entry.value < max_separation:
            index_vsx = idx[index_star_catalog]
            # if it's a new vsx star or a better match than the last match, write into dict
            if (
                index_vsx not in results_dict
                or results_dict[index_vsx].sep > entry.value
            ):
                star_id = star_descriptions[index_star_catalog].local_id
                results_dict[index_vsx] = VsxInfo(star_id, index_vsx, entry.value)
                logging.debug(f"Adding {results_dict[index_vsx]} to VSX results")
    cachedict = utils.get_localid_to_sd_dict(star_descriptions)
    # loop over dict and add the new vsx matches to the star descriptions
    for keys, vsxinfo in results_dict.items():
        logging.debug(
            f"len sd is {len(star_descriptions)}, vsxinfo.starid is {vsxinfo.starid_0}"
        )
        _add_vsx_metadata_to_star_description(
            "VSX", cachedict[vsxinfo.starid_0], vsx_dict, vsxinfo.vsx_id, vsxinfo.sep
        )
        result_ids.append(vsxinfo.starid_0)
    logging.debug(f"Added {len(results_dict)} vsx stars.")
    return star_descriptions, result_ids


# star_catalog = create_star_descriptions_catalog(star_descriptions)
def get_starid_1_for_radec(ra_deg, dec_deg, all_star_catalog, max_separation=0.01):
    star_catalog = create_generic_astropy_catalog(ra_deg, dec_deg)
    idx, d2d, d3d = match_coordinates_sky(star_catalog, all_star_catalog)
    logging.debug(
        f"get_starid_1_for_radec: len(idx):{len(idx)}, min(idx): {np.min(idx)}, max(idx): {np.max(idx)}, min(d2d): {np.min(d2d)}, "
        f"max(d2d): {np.max(d2d)}, star_id: {idx[0] + 1}, index in all_star_catalog: {all_star_catalog[idx[0]]}"
    )
    return idx[0] + 1


################ CATALOG related functions #########################


def add_metadata_to_star_descriptions(
    stars: List[StarDescription], metadata: List[object], strict: True
) -> List[StarDescription]:
    """
    Adds metadata to a list of star descriptions for later filtering
    :param strict:
    :param stars:
    :param metadata: list of strings and StarMetaData objects
    :return: the list of stars, augmented with metadata
    """
    if metadata is None:
        raise ValueError("Fill in the metadata list please.")
    for sd in stars:
        for entry in metadata:
            if isinstance(entry, str):
                sd.set_metadata(StarMetaData(key=entry), strict)
            else:
                sd.set_metadata(entry, strict)
    return stars


def _add_vsx_metadata_to_star_description(
    catalog_name: str, star: StarDescription, vsx_dict, index_vsx, separation
):
    assert star.metadata is not None
    vsx_name = vsx_dict["extradata"][index_vsx]["Name"]
    star.aavso_id = vsx_name
    match = CatalogData(
        key=catalog_name,
        catalog_id=vsx_name,
        name=vsx_name,
        separation=separation,
        coords=SkyCoord(
            vsx_dict["ra_deg_np"][index_vsx],
            vsx_dict["dec_deg_np"][index_vsx],
            unit="deg",
        ),
        extradata=vsx_dict["extradata"][index_vsx],
    )
    star.metadata = match
    return match


################ CATALOG related functions #########################

# Takes in a list of known variables and maps them to the munipack-generated star numbers
# usage:
# vsx = getVSX(settings.basedir+'SearchResults.csv')
# detections = reading.read_world_positions(settings.worldposdir)
# returns { 'name of VSX variable': [VSX_var_SkyCoord, best_starfit, best_separation] }
def find_star_for_known_vsx(vsx, detections_catalog, max_separation=0.01):
    result = {}
    logging.info(f"Searching best matches with max separation: {max_separation} ...")
    kdtree_cache = ""
    for variable in vsx:
        logging.info(f"Searching for {variable[0]}")
        idx, d2d, d3d = match_coordinates_sky(
            variable[1], detections_catalog, storedkdtree=kdtree_cache
        )
        if d2d.degree < max_separation:
            result[variable[0]] = [variable[1], idx + 1, d2d.degree]
            logging.info(f"Found result for {variable[0]} : {result[variable[0]]}")
    logging.info(f"Found {len(result)} matches")
    return result


def create_generic_astropy_catalog(ra_deg_np, dec_deg_np):
    """ Creates a SkyCoord catalog with an array of ra's and dec's. Works also with single ra dec """
    logging.debug(
        "Creating generic astropy Catalog with " + str(len(ra_deg_np)) + " objects..."
    )
    return SkyCoord(ra=ra_deg_np, dec=dec_deg_np, unit="deg")


def create_vsx_astropy_catalog(vsx_catalog_location):
    vsx_dict = vsx_pickle.read(vsx_catalog_location)
    logging.info(
        f"Creating VSX star catalog with {len(vsx_dict['ra_deg_np'])} stars using '{vsx_catalog_location}'"
    )
    return (
        create_generic_astropy_catalog(vsx_dict["ra_deg_np"], vsx_dict["dec_deg_np"]),
        vsx_dict,
    )


def create_star_descriptions_catalog(star_descriptions):
    logging.debug(
        "Creating star_descriptions star catalog with {} stars...".format(
            len(star_descriptions)
        )
    )
    ra2 = np.array([])
    dec2 = np.array([])
    for entry in star_descriptions:
        ra2 = np.append(ra2, [entry.coords.ra.deg])
        dec2 = np.append(dec2, [entry.coords.dec.deg])
    return create_generic_astropy_catalog(ra2, dec2)


# take a star_description and add some info to it: vmag, error vmag, catalog information
def add_catalog_data_to_sd(
    star, vmag, vmag_err, catalog_id, catalog_name, coord_catalog
):
    star_to_catalog_dist = star.coords.separation(coord_catalog)
    star.metadata = CatalogData(
        key=catalog_name,
        catalog_id=catalog_id,
        name=catalog_id,
        coords=coord_catalog,
        separation=star_to_catalog_dist,
        vmag=vmag,
        vmag_err=vmag_err,
    )
    logging.debug(
        "Add info: Star {} has vmag={}, error={}, dist={}".format(
            star.local_id, star.vmag, star.vmag_err, star_to_catalog_dist
        )
    )
    return star

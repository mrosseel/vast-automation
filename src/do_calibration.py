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
from astroquery.vizier import Vizier
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
    data = hdulist[0].data.astype(float)
    header = hdulist[0].header
    wcs = WCS(header)
    return wcs


def calibrate():
    w = WCS(settings.reference_header)
    return w


def old_find_reference_frame_index():
    the_dir = os.listdir(init.reference_dir)
    the_dir.sort()
    reference_frame_index = the_dir.index(init.reference_frame)
    assert the_dir[reference_frame_index] == init.reference_frame
    return reference_frame_index


# Either reads the previously chosen reference frame, or determines one if it's missing
def get_reference_frame(file_limit, reference_method):
    # Calculate or retrieve reference frame and its index
    try:
        reference_file, reference_frame_fits, reference_frame_index = reading.read_reference_frame()
        reference_frame = utils.find_file_for_index(settings.convfitsdir, reference_frame_index, '*.fts')
        logging.info(
            f"Loaded existing reference file '{reference_file}, found fits reference: {reference_frame_fits} == converted fits:{reference_frame}")
    except:
        reference_frame = reference_method(file_limit)
        reference_frame_index = utils.find_index_of_file(settings.convfitsdir, reference_frame, '*.fts')
        lines = [reference_frame, str(reference_frame_index)]
        with open(settings.basedir + 'reference_frame.txt', 'w') as f:
            f.write('\n'.join(lines))
        logging.info(
            f"Loaded reference file: '{reference_file}', found reference frame: {reference_frame}, index: {reference_frame_index}")
    assert reference_frame_index == utils.find_index_of_file(settings.convfitsdir, reference_frame, '*.fts')
    return [reference_frame, settings.convfitsdir + reference_frame, reference_frame_index]


# returns the converted_fits file with the highest compressed filesize, limited to 'limit'
# if limit is higher than listdir length, no harm
def select_reference_frame_gzip(limit):
    import gzip
    count = 0
    result = {}
    for filename in os.listdir(settings.convfitsdir):
        if count < limit:
            count = count + 1
            with open(settings.convfitsdir + filename, 'rb') as f_in:
                length = len(gzip.compress(f_in.read()))
                result[filename] = length
                logging.info(f"select reference frame length: {length}")

    sorted_by_value = sorted(result.items(), key=lambda kv: kv[1], reverse=True)
    return sorted_by_value[0][0]


# not used atm
def find_target_star(target_ra_deg, target_dec_deg, nr_results):
    target = SkyCoord(target_ra_deg, target_dec_deg, unit='deg')
    result_dict = reading.read_world_positions(settings.worldposdir)
    distances_dict = {}
    for key in result_dict:
        distances_dict[key] = target.separation(SkyCoord(result_dict[key][0], result_dict[key][1], unit='deg')).degree
    df = pd.DataFrame(list(distances_dict.items()), columns=['star_nr', 'deg_separation'])
    df.sort_values(by='deg_separation', inplace=True)
    return df[:nr_results]


############# star description utils #################

# returns list of star descriptions
def get_star_descriptions(star_id_list=None):
    # returns {'name': [ra.deg, dec.deg ]}
    positions = reading.read_world_positions(settings.worldposdir)
    result = []
    plist = "all stars"
    if star_id_list is not None:
        plist = star_id_list

    logging.info(f'Reading star descriptions for: {plist} with size {len(init.star_list)}')
    for key in positions:
        star_id = reading.filename_to_star(str(key))
        if star_id_list is None or star_id in star_id_list:
            result.append(StarDescription(local_id=reading.filename_to_star(str(key)),
                                          coords=SkyCoord(positions[key][0], positions[key][1], unit='deg')))
    return result


def select_star_descriptions(star_id_list: List[int], stars: List[StarDescription]):
    return [x for x in stars if x.local_id in star_id_list]


# returns list of star descriptions
def get_empty_star_descriptions(star_id_list=None):
    # returns {'name': [ra.deg, dec.deg ]}
    result = []
    for star in star_id_list:
        result.append(StarDescription(local_id=star))
    return result


# for testinsg
def get_random_star_descriptions(nr=10):
    result = []
    for idx in range(nr):
        result.append(StarDescription(local_id=idx,
                                      coords=SkyCoord(idx, idx, unit='deg')))
    return result


# given a list of SD's and a list of star id's, output the matching SD's in log.debug
def log_star_descriptions(star_descriptions: StarDescription, star_id_list):
    logging.debug(f"Logging {len(star_id_list)} star  descriptions:")
    for star_id in star_id_list:
        logging.debug(star_descriptions[star_id])


############# star description utils #################


# returns new list of StarDescription s with filled in local_id, upsilon match, coord
def get_candidates(threshold_prob=0.5, check_flag=False):
    result = []
    df = get_upsilon_candidates_raw(threshold_prob, check_flag)
    if df is None:
        return result
    positions = reading.read_world_positions(settings.worldposdir)
    for index, row in df.iterrows():
        upsilon_match = UpsilonData(metadata_id='Upsilon', var_type=row['label'], probability=row['probability'],
                                    flag=row['flag'], period=row['period'])
        result.append(
            StarDescription(local_id=index, metadata=upsilon_match,
                            coords=SkyCoord(positions[int(index)][0], positions[int(index)][1], unit='deg')))
    return result


# adds UpsilonMatch to existing StarDescription
def add_candidates_to_star_descriptions(stars: List[StarDescription], threshold_prob=0.5, check_flag=False):
    df = get_upsilon_candidates_raw(threshold_prob, check_flag)
    if df is None:
        return stars
    for index, row in df.iterrows():
        upsilon_match = UpsilonData(metadata_id='Upsilon', var_type=row['label'], probability=row['probability'],
                                    flag=row['flag'], period=row['period'])
        stars[index].metadata.append(upsilon_match)
    return stars


# common upsilon code
def get_upsilon_candidates_raw(threshold_prob, check_flag):
    upsilon_file = settings.basedir + 'upsilon_output.txt'
    try:
        df = pd.read_csv(upsilon_file, index_col=0)
    except:
        logging.error(f'could not read {upsilon_file}, skipping upsilon candidate processing')
        return None
    df.sort_values(by='probability', ascending=False)
    df = df[df['label'] != 'NonVar']
    df = df[df["probability"] > threshold_prob]
    if check_flag: df = df[df["flag"] != 1]
    return df


# returns {'star_id': [label, probability, flag, period, SkyCoord, match_name, match_skycoord,
# match_type, separation_deg]}
def add_vsx_names_to_star_descriptions(star_descriptions: List[StarDescription], vsxcatalogdir: str,
                                       max_separation=0.01):
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
    VsxInfo = namedtuple('VsxInfo', 'starid_0 vsx_id sep')
    for index_star_catalog, entry in enumerate(d2d):
        if entry.value < max_separation:
            index_vsx = idx[index_star_catalog]
            # if it's a new vsx star or a better match than the last match, write into dict
            if index_vsx not in results_dict or results_dict[index_vsx].sep > entry.value:
                star_id = star_descriptions[index_star_catalog].local_id
                results_dict[index_vsx] = VsxInfo(star_id, index_vsx, entry.value)
                logging.debug(f"Adding {results_dict[index_vsx]} to VSX results")
    cachedict = utils.get_star_description_cache(star_descriptions)
    # loop over dict and add the new vsx matches to the star descriptions
    for keys, vsxinfo in results_dict.items():
        logging.debug(f"len sd is {len(star_descriptions)}, vsxinfo.starid is {vsxinfo.starid_0}")
        _add_vsx_metadata_to_star_description('VSX', cachedict[vsxinfo.starid_0], vsx_dict,
                                              vsxinfo.vsx_id, vsxinfo.sep)
        result_ids.append(vsxinfo.starid_0)
    logging.debug(f"Added {len(results_dict)} vsx stars.")
    return star_descriptions, result_ids


# star_catalog = create_star_descriptions_catalog(star_descriptions)
def get_starid_1_for_radec(ra_deg, dec_deg, all_star_catalog, max_separation=0.01):
    star_catalog = create_generic_astropy_catalog(ra_deg, dec_deg)
    idx, d2d, d3d = match_coordinates_sky(star_catalog, all_star_catalog)
    logging.debug(
        f"get_starid_1_for_radec: len(idx):{len(idx)}, min(idx): {np.min(idx)}, max(idx): {np.max(idx)}, min(d2d): {np.min(d2d)}, "
        f"max(d2d): {np.max(d2d)}, star_id: {idx[0] + 1}, index in all_star_catalog: {all_star_catalog[idx[0]]}")
    return idx[0] + 1


# returns StarDescription array
def get_vsx_in_field(star_descriptions, max_separation=0.01):
    logging.info("Get VSX in field star descriptions")
    vsx_catalog, vsx_dict = create_vsx_astropy_catalog(settings.vsxcatalogdir)
    star_catalog = create_star_descriptions_catalog(star_descriptions)
    idx, d2d, d3d = match_coordinates_sky(vsx_catalog, star_catalog)
    result = []
    for index_vsx, entry in enumerate(d2d):
        if entry.value < max_separation:
            star_local_id = idx[index_vsx] + 1
            star_coords = star_catalog[star_local_id - 1]
            result_entry = StarDescription()
            result_entry.local_id = star_local_id
            result_entry.coords = star_coords
            _add_vsx_metadata_to_star_description('VSX', result_entry, vsx_dict, index_vsx, entry.value)
            result.append(result_entry)
    logging.info("Found {} VSX stars in field: {}".format(len(result), [star.local_id for star in result]))
    return result


################ CATALOG related functions #########################

def add_metadata_to_star_descriptions(stars: List[StarDescription],
                                      metadata: List[object] = None) -> List[StarDescription]:
    """
    Adds metadata to a list of star descriptions for later filtering
    :param stars:
    :param metadata: list of strings and StarMetaData objects
    :return: the list of stars, augmented with metadata
    """
    if metadata is None:
        metadata = ["SELECTED"]
    for sd in stars:
        for entry in metadata:
            if isinstance(entry, str):
                sd.metadata = StarMetaData(metadata_id=entry)
            else:
                sd.metadata.append(entry)
    return stars


def _add_vsx_metadata_to_star_description(catalog_name: str, star: StarDescription, vsx_dict, index_vsx,
                                          separation):
    assert star.metadata is not None
    vsx_name = vsx_dict['metadata'][index_vsx]['Name']
    star.aavso_id = vsx_name
    match = CatalogData(metadata_id=catalog_name, catalog_id=vsx_name,
                        name=vsx_name, separation=separation,
                        coords=SkyCoord(vsx_dict['ra_deg_np'][index_vsx], vsx_dict['dec_deg_np'][index_vsx],
                                        unit='deg'))
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
        idx, d2d, d3d = match_coordinates_sky(variable[1], detections_catalog, storedkdtree=kdtree_cache)
        if d2d.degree < max_separation:
            result[variable[0]] = [variable[1], idx + 1, d2d.degree]
            logging.info(f"Found result for {variable[0]} : {result[variable[0]]}")
    logging.info(f"Found {len(result)} matches")
    return result


def create_generic_astropy_catalog(ra_deg_np, dec_deg_np):
    logging.debug("Creating generic astropy Catalog with " + str(len(ra_deg_np)) + " objects...")
    return SkyCoord(ra=ra_deg_np, dec=dec_deg_np, unit='deg')


def create_vsx_astropy_catalog(vsx_catalog_location):
    logging.info("Creating vsx star catalog...")
    vsx_dict = vsx_pickle.read(vsx_catalog_location)
    return create_generic_astropy_catalog(vsx_dict['ra_deg_np'], vsx_dict['dec_deg_np']), vsx_dict


# NO LONGER USED
def create_upsilon_astropy_catalog(threshold_prob_candidates=0.5):
    logging.info("Creating upsilon star catalog...")
    ra2 = np.array([])
    dec2 = np.array([])
    candidates_array = get_candidates(threshold_prob_candidates)
    for candidate in candidates_array:
        ra2 = np.append(ra2, [candidate.coord.ra.deg])
        dec2 = np.append(dec2, [candidate.coord.dec.deg])
    return create_generic_astropy_catalog(ra2, dec2), candidates_array


def create_star_descriptions_catalog(star_descriptions):
    logging.info("Creating star_descriptions star catalog with {} stars...".format(len(star_descriptions)))
    ra2 = np.array([])
    dec2 = np.array([])
    for entry in star_descriptions:
        ra2 = np.append(ra2, [entry.coords.ra.deg])
        dec2 = np.append(dec2, [entry.coords.dec.deg])
    return create_generic_astropy_catalog(ra2, dec2)


def add_apass_to_star_descriptions(star_descriptions, radius=0.01, row_limit=2):
    logging.info(f"apass input {len(star_descriptions)}")
    radius_angle = Angle(radius, unit=u.deg)
    for star in star_descriptions:
        apass = get_apass_field(star.coords, radius=radius_angle, row_limit=row_limit)
        if apass is None:
            logging.info("More/less results received from APASS than expected: {}".format(
                apass.shape[0] if not apass is None and not apass.shape is None else 0))
            logging.info(apass)
            continue
        else:
            star.vmag = apass['Vmag'][0]
            star.e_vmag = apass['e_Vmag'][0]
            catalog_id = 'TODO'  # apass['recno'][0].decode("utf-8")
            coord_catalog = SkyCoord(apass['RAJ2000'], apass['DEJ2000'], unit='deg')
            mindist = star.coords.separation(SkyCoord(apass['RAJ2000'], apass['DEJ2000'], unit='deg'))
            star.metadata = CatalogData(metadata_id="APASS", catalog_id=catalog_id,
                                        name=catalog_id, coords=coord_catalog, separation=mindist)
            logging.info(
                "APASS: Star {} has vmag={}, error={:.5f}, dist={}".format(star.local_id, star.vmag, star.e_vmag,
                                                                           mindist))
    return star_descriptions


def add_vizier_ucac4_to_star_descriptions(star_descriptions: List[StarDescription], nr_threads: int, radius=0.01):
    logging.info("Retrieving ucac4 for {} star(s)".format(len(star_descriptions)))
    radius_angle = Angle(radius, unit=u.deg)

    pool = Pool(nr_threads * 2)
    func = partial(add_vizier_ucac4_to_star, radius_angle=radius_angle)
    result = []
    for entry in pool.map(func, star_descriptions):
        result.append(entry)
    return result


def add_vizier_ucac4_to_star(star: StarDescription, radius_angle):
    vizier_results = get_ucac4_field(star.coords, radius=radius_angle, row_limit=1)
    if vizier_results is None:
        logging.warning("No vizier results for star", star.local_id)
        return star
    else:
        # print('vizier results', vizier_results)
        logging.debug(f"Adding vmag: {vizier_results['Vmag'].iloc(0)[0]}")
        coord_catalog = SkyCoord(vizier_results['RAJ2000'], vizier_results['DEJ2000'], unit='deg')
        add_info_to_star_description(star, vizier_results['Vmag'].iloc(0)[0], vizier_results['e_Vmag'][0],
                                     vizier_results['UCAC4'][0].decode("utf-8"), "UCAC4", coord_catalog)
        # print(vizier_results.describe())
        # print(vizier_results.info())
        # import code
        # code.InteractiveConsole(locals=dict(globals(), **locals())).interact()
    return star


# take a star_description and add some info to it: vmag, error vmag, catalog information
def add_info_to_star_description(star, vmag, e_vmag, catalog_id, catalog_name, coord_catalog):
    star.vmag = vmag
    star.e_vmag = e_vmag
    mindist = star.coords.separation(coord_catalog)
    star.metadata = CatalogData(metadata_id=catalog_name, catalog_id=catalog_id,
                                name=catalog_id, coords=coord_catalog, separation=mindist)
    logging.debug("Add info: Star {} has vmag={}, error={}, dist={}".format(star.local_id, star.vmag, star.e_vmag,
                                                                            mindist))
    return star


def get_apass_row_to_star_descriptions(row):
    return StarDescription(coords=SkyCoord(row['RAJ2000'], row['DEJ2000'], unit='deg'),
                           vmag=row['Vmag'], e_vmag=row['e_Vmag'])


def get_ucac4_row_to_star_descriptions(row):
    ucac4 = get_apass_row_to_star_descriptions(row)
    ucac4.aavso_id = row['UCAC4']
    return ucac4


def get_apass_star_descriptions(center_coord, radius, row_limit=2):
    return get_vizier_star_descriptions(get_apass_field, get_apass_row_to_star_descriptions, center_coord,
                                        radius, row_limit)


def get_ucac4_star_descriptions(center_coord, radius, row_limit=2):
    return get_vizier_star_descriptions(get_ucac4_field, get_ucac4_row_to_star_descriptions, center_coord,
                                        radius, row_limit)


def get_vizier_star_descriptions(field_method, to_star_descr_method, center_coord, radius, row_limit=2):
    field = field_method(center_coord, radius, row_limit)
    result = []
    for index, row in field.iterrows():
        result.append(to_star_descr_method(row))
    return result


def get_apass_field(center_coord, radius, row_limit=2):
    return get_vizier_field(center_coord, radius, 'II/336', row_limit)


def get_ucac4_field(center_coord, radius, row_limit=2):
    return get_vizier_field(center_coord, radius, 'I/322A', row_limit)


# object_name = 'UCAC4 231-154752'
def get_ucac4_id_as_dataframe(object_name):
    result = Vizier.query_object(object_name, catalog=['I/322A'])
    if len(result) == 0:
        logging.warning(f"Returned result from vizier for object {object_name} is empty.")
    df = result[0].to_pandas()
    df2 = df.loc[df['UCAC4'] == object_name.split()[1].encode('UTF-8')]
    return df2


def get_vizier_field(center_coord, radius, catalog, row_limit=2):
    # Select all columns (this is a generic function so has to be all)
    v = Vizier(columns=['all'])

    # Limit the number of rows returned to row_limit
    v.ROW_LIMIT = row_limit

    result = v.query_region(center_coord,
                            radius=radius,
                            catalog=[catalog])
    try:
        df = result[0].to_pandas() if len(result) > 0 else None
    except:
        logging.error("Error processing vizier results in do_calibration")
        return None
    return df


if __name__ == '__main__':
    star_descriptions = do_calibration.get_star_descriptions(init.star_list)
    star_catalog = create_star_descriptions_catalog(star_descriptions)
    result = get_starid_1_for_radec(star_catalog, [271.234735], [-43.845816])
    logging.info(f"calibration result: {result}")

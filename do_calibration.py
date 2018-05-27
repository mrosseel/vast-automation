import os
import init
import reading
from star_description import StarDescription
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy import units as u
from astropy.coordinates import Angle
from astroquery.vizier import Vizier
import pandas as pd
import numpy as np
import vsx_pickle
import copy

def get_wcs(wcs_file):
    hdulist = fits.open(wcs_file)
    data = hdulist[0].data.astype(float)
    header = hdulist[0].header
    wcs = WCS(header)
    return wcs

def calibrate():
    w = WCS(init.reference_header)
    # xpos, ypos = w.all_world2pix(init.ra_deg, init.dec_deg, 0, ra_dec_order=True)
    # print(xpos, ypos)
    # print(w)
    return w

def find_reference_frame_index():
    the_dir = os.listdir(init.reference_dir)
    the_dir.sort()
    reference_frame_index = the_dir.index(init.reference_frame)
    assert the_dir[reference_frame_index] == init.reference_frame
    return reference_frame_index


# returns 'path + phot????.pht', the photometry file matched with the reference frame
def find_reference_photometry(reference_frame_index):
    the_dir = os.listdir(init.photometrydir)
    the_dir.sort()
    return init.photometrydir + the_dir[reference_frame_index]


def find_reference_matched(reference_frame_index):
    the_dir = os.listdir(init.matchedphotometrydir)
    the_dir.sort()
    return init.matchedphotometrydir + the_dir[reference_frame_index]


def find_target_star(target_ra_deg, target_dec_deg, nr_results):
    target = SkyCoord(target_ra_deg, target_dec_deg, unit='deg')
    result_dict = reading.read_world_positions(init.worldposdir)
    distances_dict = {}
    for key in result_dict:
        distances_dict[key] = target.separation(SkyCoord(result_dict[key][0], result_dict[key][1], unit='deg')).degree
    df = pd.DataFrame(list(distances_dict.items()), columns=['star_nr', 'deg_separation'])
    df.sort_values(by='deg_separation', inplace=True)
    return df[:nr_results]


def find_target_stars(max_deg_separation):
    target = SkyCoord(target_ra_deg, target_dec_deg, unit='deg')
    result_dict = reading.read_world_positions(init.worldposdir)
    distances_dict = {}
    for key in result_dict:
        distances_dict[key] = target.separation(SkyCoord(result_dict[key][0], result_dict[key][1], unit='deg')).degree
    df = pd.DataFrame(list(distances_dict.items()), columns=['filename', 'deg_separation'])
    df.sort_values(by='deg_separation', inplace=True)
    return df[:nr_results]


# returns StarDescription with filled in local_id, upsilon, coord
def get_candidates(threshold_prob=0.5, check_flag=False):
    df = pd.DataFrame.from_csv(init.basedir + 'upsilon_output.txt')
    df.sort_values(by='probability', ascending=False)
    df = df[df['label'] != 'NonVar']
    df = df[df["probability"] > threshold_prob]
    if check_flag: df = df[df["flag"] != 1]
    positions = reading.read_world_positions(init.worldposdir)
    result = []
    for index, row in df.iterrows():
        # print(index, row)
        result.append(
            StarDescription(local_id=index, upsilon=(row['label'], row['probability'], row['flag'], row['period']),
                            coords=SkyCoord(positions[int(index)][0], positions[int(index)][1], unit='deg')))
    return result

# returns StarDescription with filled in local_id, upsilon, coord
def add_upsilon_data(star_descriptions):
    df = pd.DataFrame.from_csv(init.basedir + 'upsilon_output.txt')
    for star in star_descriptions:
        try:
            row = df.iloc[df.index.get_loc(star.local_id)]
            star.upsilon = (row['label'], row['probability'], row['flag'], row['period'])
        except:
            continue
    return star_descriptions

# returns list of star descriptions
def get_star_descriptions(starlist=None):
    # returns {'name': [ra.deg, dec.deg ]}
    positions = reading.read_world_positions(init.worldposdir)
    result = []
    print('Reading star descriptions for:', starlist if not None else 'all stars')
    for key in positions:
        star_id = reading.filename_to_star(str(key))
        if starlist is None or star_id in starlist:
            result.append(StarDescription(local_id=reading.filename_to_star(str(key)), coords=SkyCoord(positions[key][0], positions[key][1], unit='deg')))
    return result

# returns [index, skycoord, type]
def get_VSX(the_file):
    raise NotImplementedError
    result = []
    df = pd.DataFrame.from_csv(the_file)
    # print(df.head())
    for index, row in df.iterrows():
        skycoord = SkyCoord(row['Coords'], unit=(u.hourangle, u.deg))
        result.append(StarDescription(match={'index': index, 'type': row['Type']}, coord=skycoord))
    return result


# returns {'star_id': [label, probability, flag, period, SkyCoord, match_name, match_skycoord, match_type, separation_deg]}
def add_vsx_names_to_star_descriptions(star_descriptions, threshold_prob_candidates=0.5, max_separation=0.01):
    print("Adding VSX names to star descriptions")
    result = star_descriptions  # no deep copy for now
    # copy.deepcopy(star_descriptions)
    vsx_catalog, vsx_dict = create_vsx_astropy_catalog()
    star_catalog = create_star_descriptions_catalog(star_descriptions)
    idx, d2d, d3d = match_coordinates_sky(star_catalog, vsx_catalog)
    print(len(idx))
    found = 0
    for index, entry in enumerate(d2d):
        if entry.value < max_separation:
            vsx_index = idx[index]
            result[index].match = ('VSX', entry.value,
                                   {'name': vsx_dict['metadata'][vsx_index]['Name'],
                                    'coords': SkyCoord(vsx_dict['ra_deg_np'][index], vsx_dict['dec_deg_np'][vsx_index],
                                                       unit='deg')})
            result[index].aavso_id = vsx_dict['metadata'][vsx_index]['Name']
            found = found + 1
    return result


# returns StarDescription array
def get_vsx_in_field(star_descriptions, max_separation=0.01):
    print("Adding VSX names to star descriptions")
    vsx_catalog, vsx_dict = create_vsx_astropy_catalog()
    star_catalog = create_star_descriptions_catalog(star_descriptions)
    idx, d2d, d3d = match_coordinates_sky(vsx_catalog, star_catalog)
    result = []
    for index, entry in enumerate(d2d):
        if entry.value < max_separation:
            star_local_id = idx[index]
            star_coords = star_catalog[star_local_id]
            vsx_coords = SkyCoord(vsx_dict['ra_deg_np'][index], vsx_dict['dec_deg_np'][index], unit='deg')
            result_entry = StarDescription()
            result_entry.local_id = star_local_id
            result_entry.coords = star_coords
            result_entry.aavso_id = vsx_dict['metadata'][index]['Name']
            result_entry.match = ('VSX', entry.value,
                                  {'name': result_entry.aavso_id,
                                   'coords': vsx_coords})
            result.append(result_entry)
    print("Found {} stars".format(len(result)))
    return result

# count the number of vsx in the field
def count_vsx_in_field():
    #bla
    print('bla')


# Takes in a list of known variables and maps them to the munipack-generated star numbers
# usage:
# vsx = getVSX(init.basedir+'SearchResults.csv')
# detections = reading.read_world_positions(init.worldposdir)
# returns { 'name of VSX variable': [VSX_var_SkyCoord, best_starfit, best_separation] }
def find_star_for_known_vsx(vsx, detections_catalog, max_separation=0.01):
    result = {}
    print("Searching best matches with max separation", max_separation, "...")
    kdtree_cache = ""
    for variable in vsx:
        print("Searching for", variable[0])
        idx, d2d, d3d = match_coordinates_sky(variable[1], detections_catalog, storedkdtree=kdtree_cache)
        if d2d.degree < max_separation:
            result[variable[0]] = [variable[1], idx + 1, d2d.degree]
            print("Found result for ", variable[0], ":", result[variable[0]])
    print("Found", len(result), "matches")
    return result


def create_generic_astropy_catalog(ra_deg_np, dec_deg_np):
    print("Creating astropy Catalog with " + str(len(ra_deg_np)) + " objects...")
    return SkyCoord(ra=ra_deg_np, dec=dec_deg_np, unit='deg')


def create_detections_astropy_catalog(detections):
    print("Created detected stars catalog...")
    ra2 = np.array([])
    dec2 = np.array([])
    # rewrite ra/dec to skycoord objects
    for key in detections:
        ra2 = np.append(ra2, [detections[key][0]])
        dec2 = np.append(dec2, [detections[key][1]])
    return create_generic_astropy_catalog(ra2, dec2)


def create_vsx_astropy_catalog():
    print("Creating vsx star catalog...")
    vsx_dict = vsx_pickle.read(init.vsx_catalog_path)
    return create_generic_astropy_catalog(vsx_dict['ra_deg_np'], vsx_dict['dec_deg_np']), vsx_dict


def create_upsilon_astropy_catalog(threshold_prob_candidates=0.5):
    print("Creating upsilon star catalog...")
    ra2 = np.array([])
    dec2 = np.array([])
    candidates_array = get_candidates(threshold_prob_candidates)
    for candidate in candidates_array:
        ra2 = np.append(ra2, [candidate.coord.ra.deg])
        dec2 = np.append(dec2, [candidate.coord.dec.deg])
    return create_generic_astropy_catalog(ra2, dec2), candidates_array


def create_star_descriptions_catalog(star_descriptions):
    print("Creating star_descriptions star catalog with {} stars...".format(len(star_descriptions)))
    ra2 = np.array([])
    dec2 = np.array([])
    for entry in star_descriptions:
        ra2 = np.append(ra2, [entry.coords.ra.deg])
        dec2 = np.append(dec2, [entry.coords.dec.deg])
    return create_generic_astropy_catalog(ra2, dec2)


def add_apass_to_star_descriptions(star_descriptions, radius=0.01, row_limit=2):
    print("apass input", len(star_descriptions))
    radius_angle = Angle(radius, unit=u.deg)
    for star in star_descriptions:
        apass = get_apass_field(star.coords, radius=radius_angle, row_limit=row_limit)
        if apass is None:
            print("More/less results received from APASS than expected: {}".format(apass.shape[0] if not apass is None and not apass.shape is None else 0))
            print(apass)
            continue
        else:
            star.vmag = apass['Vmag'][0]
            star.e_vmag = apass['e_Vmag'][0]
            mindist = star.coords.separation(SkyCoord(apass['RAJ2000'], apass['DEJ2000'], unit='deg'))
        print("Star {} has vmag={}, error={:.5f}, dist={}".format(star.local_id, star.vmag, star.e_vmag, mindist))
    return star_descriptions

def add_ucac4_to_star_descriptions(star_descriptions, radius=0.01):
    print("apass input", len(star_descriptions))
    radius_angle = Angle(radius, unit=u.deg)
    for star in star_descriptions:
        vizier_results = get_ucac4_field(star.coords, radius=radius_angle, row_limit=1)
        if vizier_results is None:
            print("More/less results received from APASS than expected: {}".format(vizier_results.shape[0] if not vizier_results is None and not vizier_results.shape is None else 0))
            print(vizier_results)
            continue
        else:
            print('vizier results', vizier_results, vizier_results['Vmag'][0])
            star.vmag = vizier_results['Vmag'][0]
            star.e_vmag = vizier_results['e_Vmag'][0]
            mindist = star.coords.separation(SkyCoord(vizier_results['RAJ2000'], vizier_results['DEJ2000'], unit='deg'))
        print(vizier_results.describe())
        print(vizier_results.info())
        print(star.e_vmag)
        print("Star {} has vmag={}, error={}, dist={}".format(star.local_id, star.vmag, star.e_vmag, mindist))
    return star_descriptions

def get_apass_row_to_star_descriptions(row):
    return StarDescription(coords=SkyCoord(row['RAJ2000'], row['DEJ2000'], unit='deg'),
                           vmag=row['Vmag'], e_vmag=row['e_Vmag'])

def get_ucac4_row_to_star_descriptions(row):
    ucac4 =  get_apass_row_to_star_descriptions(row)
    ucac4.aavso_id = row['UCAC4']
    return ucac4

def get_apass_star_descriptions(center_coord, radius, row_limit=2):
    return get_vizier_star_descriptions(get_apass_field, get_apass_row_to_star_descriptions,center_coord,
                                        radius, row_limit)

def get_ucac4_star_descriptions(center_coord, radius, row_limit=2):
    return get_vizier_star_descriptions(get_ucac4_field, get_ucac4_row_to_star_descriptions,center_coord,
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

def get_vizier_field(center_coord, radius, catalog, row_limit=2):
    # Select all columns
    v = Vizier(columns=['all'])

    # Limit the number of rows returned to row_limit
    v.ROW_LIMIT = row_limit

    result = v.query_region(center_coord,
                            radius=radius,
                            catalog=[catalog])
    df = result[0].to_pandas() if len(result) > 0 else None
    return df


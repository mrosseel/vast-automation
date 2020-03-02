import logging
import argparse
import numpy as np
import do_calibration
import utils
from astropy.coordinates import SkyCoord
from typing import List, Tuple, Dict
import star_description
from star_description import StarDescription
from star_metadata import CompStarData
from ucac4 import UCAC4
import math
from comparison_stars import ComparisonStars
from pathlib import PurePath
import operator
from pandas import DataFrame
import tqdm
import reading
from utils import StarDict



# receives ucac numbers, fetches ucac coords and compares them to world_position coords
def get_fixed_compstars(star_descriptions: List[StarDescription], comparison_stars: List[str]):
    logging.info(f"Using fixed compstars {comparison_stars}")
    ucac4 = UCAC4()
    star_ids_1 = []
    star_desc_result = []
    star_catalog = do_calibration.create_star_descriptions_catalog(star_descriptions)
    for ucac_id in comparison_stars:
        # getting star_id_1
        ucacsd = ucac4.get_ucac4_star_description_fromid(ucac_id)
        ra, dec = ucacsd.coords.ra, ucacsd.coords.dec
        star_id_1 = do_calibration.get_starid_1_for_radec([ra], [dec], star_catalog)
        star_ids_1.append(star_id_1)
        # adding info to star_description
        star = star_descriptions[star_id_1 - 1]
        logging.info(f"Compstar match: {ucacsd.aavso_id} with {star.local_id} ({ra}, {dec})")
        do_calibration.add_info_to_star_description(star, ucacsd.vmag, ucacsd.e_vmag,
                                                    ucacsd.aavso_id,
                                                    "UCAC4", SkyCoord(ra, dec, unit='deg'))
        star_desc_result.append(star)
        logging.debug(f"Using fixed compstar '{ucac_id}' with Vmag: '{ucacsd.vmag}' and star id: {star_id_1}")
    return [x.local_id for x in star_desc_result], star_desc_result


def get_calculated_compstars(vastdir, stardict: StarDict, ref_jd, maglimit=15, starlimit=1000, ):
    likely = _get_list_of_likely_constant_stars(vastdir)
    likely_sd: List[StarDescription] = [stardict[x] for x in likely if x in stardict]

    star_realmag = {}
    for star in likely_sd:
        comp_magdict = reading.read_magdict_for_star(vastdir, star.local_id)
        # logging.info(f"Read comp magdict for {star}: {read_comp_magdict}")
        if ref_jd in comp_magdict:
            star_realmag[star.local_id] = comp_magdict[ref_jd]

    # restrict to maglimit
    mag_list = list(filter(lambda x: x.vmag < maglimit, likely_sd))
    # only want the best stars with max possible observations
    max_obs_sorted = sorted(mag_list, key=lambda x: x.obs, reverse=True)
    max_obs = np.max([x.obs for x in likely_sd])
    max_obs_clipped = max_obs_sorted[:starlimit]
    logging.info(f"Picked {len(max_obs_clipped)} stars with last star having "
                 f"{max_obs_clipped[-2:-1][0].obs * 100 / max_obs:.2f} % of max observations ({max_obs})")


    def limit(array, max_size):
        return min(len(array), max_size)


    # stars are sorted according to their magnitude error
    min_err_stars = sorted(max_obs_clipped, key=lambda entry: star_realmag[entry.local_id][
        1] if entry.local_id in star_realmag else 1000)
    # only first 100 are kept
    min_err_stars_clipped = min_err_stars[:limit(min_err_stars, 100)]
    return [x.local_id for x in min_err_stars_clipped], min_err_stars_clipped




def _get_list_of_likely_constant_stars(vastdir):
    result = []
    with open(PurePath(vastdir, 'vast_list_of_likely_constant_stars.log'), 'r') as infile:
        for line in infile:
            if line is not None:
                result.append(utils.get_starid_from_outfile(line.strip()))
    return result


def calculate_ensemble_photometry(df: DataFrame, comp_stars: ComparisonStars, ensemble_method):
    assert comp_stars is not None
    logging.debug(f"Start calculate_real with {df.shape[0]} rows and {len(comp_stars.observations)} comp stars.")

    # error = sqrt((vsig**2+(1/n sum(sigi)**2)))
    realV = []
    realErr = []
    for index, row in df.iterrows():
        # logging.debug(f"row is {row}, comp_mags is {comp_stars.comp_catalogmags}")
        comp_obs = []
        comp_err = []
        comp_real = []
        for idx, compstar in enumerate(comp_stars.observations):
            jd = row['JD']
            if jd in compstar:
                comp_obs.append(compstar[jd][0])
                comp_err.append(compstar[jd][1])
                comp_real.append(comp_stars.comp_catalogmags[idx])
            else:
                logging.debug(f"Key error for {row['JD']}, {comp_stars.ids}")

        if len(comp_obs) > 0 and len(comp_obs) == len(comp_err):
            v = ensemble_method(row['Vrel'], comp_obs, comp_err, comp_real)
            err = math.sqrt(math.pow(row['err'], 2) + math.pow(np.mean(comp_err), 2))
            realV.append(v)
            realErr.append(err)
        else:  # error in the comparison stars
            logging.error(f"len comp obs: {len(comp_obs)}, len(comperr): {len(comp_err)}, obs:{comp_obs}, err:{comp_err}, index:{index} , row:{row}, df: {df.describe()}, df: {df.info()}, comp_stars: {comp_stars}")
            #realV.append(row['Vrel'])
            #realErr.append(row['err'])
    return realV, realErr


# DF should have JD, Vrel, err
def mean_value_ensemble_method(vrel, comp_obs, comp_err, comp_real):
    # Vobs - Cobs + Creal = V
    return np.mean(np.add(np.subtract(np.repeat(vrel, len(comp_obs)), comp_obs), comp_real))


def weighted_value_ensemble_method(vrel, comp_obs, comp_err, comp_real):
    vx = np.add(np.subtract(np.repeat(vrel, len(comp_obs)), comp_obs), comp_real)
    sum_vx_divided_by_errors = np.sum(np.divide(vx, comp_err))
    sum_inverse_errors = np.sum(np.divide(np.ones(len(comp_err)), comp_err))
    logging.debug(f"summags: {sum_vx_divided_by_errors}, sum inverse: {sum_inverse_errors}")
    vw = np.divide(sum_vx_divided_by_errors, sum_inverse_errors)  # eq 6
    return vw


def add_closest_compstars(stars: List[StarDescription], comp_stars: ComparisonStars, limit=10):
    assert limit > 1  # we need one extra K star for aavso ensemble output
    for star in tqdm.tqdm(stars, total=len(stars), desc='Adding closest comparison stars'):
        closest_ids = _closest_compstar_ids(star, comp_stars, limit)
        star.metadata = CompStarData(compstar_ids=closest_ids[:-1], extra_id=closest_ids[-1])


def _closest_compstar_ids(star: StarDescription, comp_stars: ComparisonStars, limit=10) -> List[int]:
    sorted_closest = sorted(comp_stars.star_descriptions, key=lambda x: x.coords.separation(star.coords))
    sorted_clipped = sorted_closest[:min(len(sorted_closest), limit)]
    sorted_ids = [x.local_id for x in sorted_clipped]
    return sorted_ids


def filter_comparison_stars(star: StarDescription, comp_stars: ComparisonStars) -> Tuple[ComparisonStars,
                                                                                         ComparisonStars]:
    compstar_match: star_description.CompStarData = star.get_metadata("COMPSTARS")
    sd_ids = compstar_match.compstar_ids
    # skip part of the work if the list is equal. Biggest set first otherwise diff is always empty
    if len(set(comp_stars.ids).difference(set(sd_ids))) == 0:
        return comp_stars
    filtered_compstars = comp_stars.get_filtered_comparison_stars(sd_ids)
    extra_compstar = comp_stars.get_filtered_comparison_stars([compstar_match.extra_id])
    logging.debug(f"get star compstars from catalog: {len(filtered_compstars.ids)}, {filtered_compstars.ids}")
    return filtered_compstars, extra_compstar

# if __name__ == '__main__':
#     logging.getLogger().setLevel(logging.DEBUG)
#     logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
#     parser = argparse.ArgumentParser(description='Getting the compstars')
#     parser.add_argument('-f', '--fixed', help="Use the comparison stars as specified in init.py", action="store_true")
#     parser.add_argument('-c', '--calculated', help="Calculate the comparison stars to use", action="store_true")
#     parser.add_argument('-d', '--descriptions', help="Calculate star descriptions", action="store_true")
#     args = parser.parse_args()
#     if args.fixed:
#         result = get_fixed_compstars()
#         logging.info(f"result: {result}")
#     elif args.calculated:
#         # stddevs, _, apertures, apertureidx, _, _, counts = do_aperture.main(the_dir=settings.matchedphotometrydir,
#         #                                                                     percentage=init.aperture_find_percentage)
#         result = get_calculated_compstars(apertureidx, stddevs, counts)
#         logging.info(f"result: {result}")
#     # if args.descriptions:
#     #     get_comparison_star_descriptions(result)

import logging
import argparse
import numpy as np
import do_calibration
import utils
from astropy.coordinates import SkyCoord
from photometry_blob import PhotometryBlob
from typing import List, Tuple, Dict
import star_description
from star_description import StarDescription
from ucac4 import UCAC4
import math
from comparison_stars import ComparisonStars
from pathlib import PurePath
import operator
import tqdm


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


def get_calculated_compstars(vastdir, stardict: Dict[int, StarDescription]):
    likely = _get_list_of_likely_constant_stars(vastdir)
    stars = [stardict[x] for x in likely if x in stardict]

    # only want the best stars
    max_obs = np.max([x.obs for x in stars])
    max_obs_stars = [x for x in stars if x.obs == max_obs]
    logging.info(f"Stars with maximum of observations: {len(max_obs_stars)}")


    def limit(array, max_size):
        return min(len(array), max_size)


    min_err_stars = sorted(max_obs_stars, key=operator.attrgetter('e_vmag'))
    min_err_stars_clipped = min_err_stars[:limit(min_err_stars, 100)]
    logging.info(f"Using {len(min_err_stars_clipped)} calculated comparison stars: {min_err_stars_clipped[:3]}")
    return [x.local_id for x in min_err_stars_clipped], min_err_stars_clipped


def _get_list_of_likely_constant_stars(vastdir):
    result = []
    with open(PurePath(vastdir, 'vast_list_of_likely_constant_stars.log'), 'r') as infile:
        for line in infile:
            if line is not None:
                result.append(utils.get_starid_from_outfile(line.strip()))
    return result


# DF should have JD, Vrel, err
def calculate_mean_value_ensemble_photometry(df, comp_stars: ComparisonStars):
    assert comp_stars is not None
    logging.debug(f"Start calculate_real with {df.shape[0]} rows and {len(comp_stars.observations)} comp stars.")

    # the average of the real magnitude of the comparison stars
    meanreal = np.mean(comp_stars.comp_catalogmags)
    logging.debug(f"meanreal is {meanreal}")
    # error = sqrt((vsig**2+(1/n sum(sigi)**2)))
    realV = []
    realErr = []
    for index, row in df.iterrows():
        # logging.debug(f"row is {row}, comp_mags is {comp_stars.comp_catalogmags}")
        comp_obs = []
        comp_err = []
        for compstar in comp_stars.observations:
            jd = row['JD']
            if jd in compstar:
                comp_obs.append(float(compstar[jd][0]))
                comp_err.append(float(compstar[jd][1]))
            else:
                logging.debug(f"Key error for {row['JD']}, {comp_stars.ids}")

        meanobs = -1
        meanerr = -1
        if len(comp_obs) > 0 and len(comp_obs) == len(comp_err):
            meanobs = np.mean(comp_obs)
            # Vobs - Cobs + Creal = V
            realV.append(row['Vrel'] - meanobs + meanreal)
            meanerr = math.sqrt(math.pow(row['err'], 2) + math.pow(np.mean(comp_err), 2))
            realErr.append(meanerr)
        else:  # error in the comparison stars
            realV.append(row['Vrel'])
            realErr.append(row['err'])
        if meanobs == -1 or meanerr == -1:
            logging.info(f"vrel: {row['Vrel']}, meanobs: {meanobs}, compobs: {comp_obs},  meanreal: {meanreal}, "
                         f"real {comp_stars.comp_catalogmags},  vrel: {row['Vrel'] - meanobs + meanreal}, meanerr: {meanerr},"
                         f"nr of compstar observations={len(compstar)}, nr of variable observations={len(df)}")
            logging.info(f"{len(comp_obs)}, {len(comp_obs)} == {len(comp_err)}")
        # logging.debug(f"Results of photometry: V diff: {df['Vrel'].mean()-np.mean(realV)}, err diff: {df['err'].mean()-np.mean(realErr)}")
    return realV, realErr


def add_closest_compstars(stars: List[StarDescription], comp_stars: ComparisonStars):
    for star in stars:
        star_description.add_compstar_match(star, closest_compstar_ids(star, comp_stars))


def closest_compstar_ids(star: StarDescription, comp_stars: ComparisonStars, limit=10) -> List[int]:
    sorted_closest = sorted(comp_stars.star_descriptions, key=lambda x: x.coords.separation(star.coords))
    sorted_clipped = sorted_closest[:min(len(sorted_closest), limit)]
    sorted_ids = [x.local_id for x in sorted_clipped]
    return sorted_ids
    # return comp_stars.get_filtered_comparison_stars(sorted_ids)


def get_star_compstars_from_catalog(star: StarDescription, comp_stars: ComparisonStars):
    compstar_match: star_description.CompStarMatch = star.get_catalog("COMPSTARS")
    sd_ids = compstar_match.compstar_ids
    # skip part of the work if the list is equal
    if len(set(sd_ids).difference(comp_stars.ids)) == 0:
        return comp_stars
    filtered_compstars = comp_stars.get_filtered_comparison_stars(sd_ids)
    return filtered_compstars


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description='Getting the compstars')
    parser.add_argument('-f', '--fixed', help="Use the comparison stars as specified in init.py", action="store_true")
    parser.add_argument('-c', '--calculated', help="Calculate the comparison stars to use", action="store_true")
    parser.add_argument('-d', '--descriptions', help="Calculate star descriptions", action="store_true")
    args = parser.parse_args()
    if args.fixed:
        result = get_fixed_compstars()
        logging.info(f"result: {result}")
    elif args.calculated:
        # stddevs, _, apertures, apertureidx, _, _, counts = do_aperture.main(the_dir=settings.matchedphotometrydir,
        #                                                                     percentage=init.aperture_find_percentage)
        result = get_calculated_compstars(apertureidx, stddevs, counts)
        logging.info(f"result: {result}")
    # if args.descriptions:
    #     get_comparison_star_descriptions(result)

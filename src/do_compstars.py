import do_calibration
import logging
import argparse
import numpy as np
from astropy.coordinates import SkyCoord
from photometry_blob import PhotometryBlob
from typing import List, Tuple
from star_description import StarDescription
import ucac4
import math
from comparison_stars import ComparisonStars


# receives ucac numbers, fetches ucac coords and compares them to world_position coords
def get_fixed_compstars(star_descriptions: List[StarDescription], comparison_stars: List[str]):
    logging.info(f"Using fixed compstars {comparison_stars}")
    star_ids_1 = []
    star_desc_result = []
    star_catalog = do_calibration.create_star_descriptions_catalog(star_descriptions)
    for ucac_id in comparison_stars:
        # getting star_id_1
        ucacsd = ucac4.get_ucac4_star_description(ucac_id)
        ra, dec = ucacsd.coords.ra, ucacsd.coords.dec
        star_id_1 = do_calibration.get_starid_1_for_radec([ra], [dec], star_catalog)
        star_ids_1.append(star_id_1)
        # adding info to star_description
        star = star_descriptions[star_id_1 - 1]
        logging.info(f"Compstar match: {ucacsd.aavso_id} with {star.local_id}")
        do_calibration.add_info_to_star_description(star, ucacsd.vmag, ucacsd.e_vmag,
                                                    ucacsd.aavso_id,
                                                    "UCAC4", SkyCoord(ra, dec, unit='deg'))
        star_desc_result.append(star)
        logging.debug(f"Using fixed compstar '{ucac_id}' with Vmag: '{ucacsd.vmag}' and star id: {star_id_1}")
    return [x.local_id for x in star_desc_result], star_desc_result


def get_calculated_compstars(apertureidx, photometry_blob: PhotometryBlob):
    stddevs = photometry_blob.stddevs
    counts = photometry_blob.counts
    selectedcounts = counts[apertureidx]

    # winners = np.argwhere(np.amax(selectedcounts)).flatten()
    logging.info(f"Before selection on counts: {len(selectedcounts)}")
    logging.info(f"counts threshold {np.max(selectedcounts) * 0.9}")
    count_mask = selectedcounts < np.max(selectedcounts) * 0.9
    logging.info(f"Count mask: {count_mask}")
    logging.info(f"After selection on counts: {np.sum(count_mask)}")
    #    logging.info(f"Counts: {sorted(selectedcounts)}")

    masked_std = np.ma.masked_array(stddevs[apertureidx], count_mask)
    logging.info(f"Masked std: {masked_std}")
    logging.info(f"masked std len: {masked_std.count()}")

    nrtopstars = min(10, len(masked_std))
    compstars_0 = masked_std.argsort(fill_value=99999)[:nrtopstars]
    logging.info(f"compstars: {compstars_0}, len: {len(compstars_0)}")

    # compstar_0 = np.argmin(stddevs[apertureidx], axis=0)
    logging.info(f"Compstars_0 with minimum stdev in the chosen aperture: {compstars_0}")
    logging.info(f"Compstars stddev: {stddevs[apertureidx][compstars_0]}")
    logging.info(f"Compstars counts: {selectedcounts[compstars_0]}")
    # for star in compstars_0:
    #     logging.info(f"Error for compstar_0:{star} is \t {errors[star].median()}")
    comparison_stars_1 = (compstars_0 + 1).tolist()
    comparison_stars_1_desc = do_calibration.get_star_descriptions(comparison_stars_1)
    return comparison_stars_1, comparison_stars_1_desc


# get star_descriptions and ucuc4 info
def get_comparison_star_descriptions(comparison_stars_1):
    logging.info(f'Getting star description of comparison star: {comparison_stars_1}')
    comp_star_description = do_calibration.get_star_descriptions(comparison_stars_1)
    comparison_stars_1_desc = do_calibration.add_ucac4_to_star_descriptions(comp_star_description)
    logging.info(f'Added apass info to comparison stars: {comparison_stars_1_desc}')

    if np.isnan(comparison_stars_1_desc[0].vmag):
        logging.info("Comparison star has nan vmag, will screw up everything coming after")
        exit()
    return comparison_stars_1_desc


def calculate_real_mag_and_err(df, comp_stars: ComparisonStars):
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
                logging.info(f"Key error for {row['JD']}, {comp_stars.ids}")

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

import do_calibration
import logging
from init_loader import init, settings
import argparse
import do_aperture
import numpy as np
from astropy.coordinates import SkyCoord
from photometry_blob import PhotometryBlob
from typing import List, Tuple
from star_description import StarDescription
import ucac4

# main entry
def select_compstars(stddevs, apertureidx, counts):
    if len(init.comparison_stars) == 0:  # if there are no given comparison stars
        compstars_1 = get_calculated_compstars(apertureidx, stddevs, counts)
    else:  # otherwise use them
        compstars_1 = get_fixed_compstars()
    return compstars_1


# reads ucac numbers from init, fetches ucac coords and compares them to world_position coords
def get_fixed_compstars(star_descriptions: List[StarDescription], comparison_stars):
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
        logging.info(f"matched {ucacsd} with {star}")
        do_calibration.add_info_to_star_description(star, ucacsd.vmag, ucacsd.e_vmag,
                                                    ucacsd.aavso_id,
                                                    "UCAC4", SkyCoord(ra, dec, unit='deg'))
        star_desc_result.append(star)
        logging.debug(f"Using fixed compstar '{ucac_id}' with Vmag: '{ucacsd.vmag}' and star id: {star_id_1}")
    return star_ids_1, star_desc_result

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
        stddevs, _, apertures, apertureidx, _, _, counts = do_aperture.main(the_dir=settings.matchedphotometrydir,
                                                                            percentage=init.aperture_find_percentage)
        result = get_calculated_compstars(apertureidx, stddevs, counts)
        logging.info(f"result: {result}")
    # if args.descriptions:
    #     get_comparison_star_descriptions(result)

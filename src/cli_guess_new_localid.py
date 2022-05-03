import argparse
import logging

from astropy.coordinates import match_coordinates_sky

import do_calibration
import utils
import utils_sd
from ucac4 import UCAC4
import main_vast

def main(vastdir, file):
    ucac4 = UCAC4()

    # construct star descriptions
    sds = utils_sd.construct_star_descriptions(vastdir, None)
    star_catalog = do_calibration.create_star_descriptions_catalog(sds)
    df = main_vast.read_localid(file)
    for idx, row in df.iterrows():
        rowucac4 = row["ucac4_name"]
        result = ucac4_to_localid(ucac4, rowucac4, star_catalog, sds)
        row["local_id"] = result 
    print(df.to_csv())
    print(results)

def ucac4_to_localid(ucac4, ucac4_name, star_catalog, sds):
    # get UCAC4 sd for search later
    ucac_details = ucac4.get_ucactuple_from_id(utils.get_full_ucac4_id(ucac4_name))
    ucacsd = ucac4.get_star_description_from_tuple(ucac_details)
    logging.info(f"Found UCAC4: {ucacsd}")

    # Get the 10 closest neighbours
    #sd_dict = utils.get_localid_to_sd_dict(sds)
    neighbours = []
    for neigh in range(1, 22):
        idx, d2d, _ = match_coordinates_sky(
            ucacsd.coords, star_catalog, nthneighbor=neigh
        )
        neighbours.append(sds[idx])
    ucac4.add_sd_metadatas(neighbours)
    # logging.info("\n" + "\n".join([f"Star {x.local_id}: {x}" for x in neighbours]))
    return neighbours[0].local_id


if __name__ == "__main__":
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(
        description="Gets a list of stars with localids and ucac names (standard CSV output of site) and generates the same file with corrected localids"
    )
    parser.add_argument(
        "-d",
        "--datadir",
        help="The directory where the data can be found (usually the vast dir)",
        nargs="?",
        required=True,
    )
    parser.add_argument("-f", "--file", help="file with localid data", required=True)
    parser.add_argument(
        "-x", "--verbose", help="Set logging to debug mode", action="store_true"
    )
    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    main(args.datadir, args.file)

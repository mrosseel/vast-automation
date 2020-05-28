from typing import List, Dict, Tuple
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from photutils import aperture_photometry, CircularAperture
from astropy.coordinates import match_coordinates_sky
from pathlib import Path
import main_vast
import utils
import os
from reading import ImageRecord
import do_calibration
import utils_sd
import subprocess
import logging
import time
import reading
import argparse
from datetime import datetime
import random
import toml
from collections import namedtuple

from star_description import StarDescription
from ucac4 import UCAC4, MinimalStarTuple

padding = 0
dpi = 100
RefFrame = namedtuple('RefFrame', 'ref_jd path_to_solved path_to_reference_frame')
ucac4 = UCAC4()


def process(stars_input):
    logging.debug(f"Stars: {stars_input}")
    stars = list(map(lambda x: utils.get_full_ucac4_id(x), stars_input))
    star1 = ucac4.get_ucactuple_from_id(stars[0])[0]
    sd1 = ucac4.get_star_description_from_tuple(*star1)
    star2 = ucac4.get_ucactuple_from_id(stars[1])[0]
    sd2 = ucac4.get_star_description_from_tuple(*star2)
    logging.info(f"Separation from {stars[0]} to {stars[1]} = {sd1.coords.separation(sd2.coords)}")


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description='munipack automation cli')
    parser.add_argument('-s', '--stars', help="List the star id's to plot", nargs='+', required=False)
    parser.add_argument('-x', '--verbose', help="Set logging to debug mode", action="store_true")
    args = parser.parse_args()
    # add the handlers to the logger
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    process(args.stars)

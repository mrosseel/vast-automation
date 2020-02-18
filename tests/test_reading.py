# from .context import src
import unittest
import sys

import star_metadata
import utils
from star_description import StarDescription
from astropy.coordinates import SkyCoord
import logging
import random
import os
from pathlib import PurePath
import logging
import pandas as pd
import numpy as np
from star_metadata import CatalogData
import reading

logging.getLogger().setLevel(logging.DEBUG)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")

test_file_path = PurePath(os.getcwd(), 'tests', 'data')


class TestReading(unittest.TestCase):

    def setUp(self) -> None:
        self.a = self.stardesc(1, 1, 1)
        self.b = self.stardesc(2, 2, 2)
        self.c = self.stardesc(3, 3, 3)
        self.a.path = './tests/data/out02391.dat'
        self.b.path = './tests/data/out02391.dat'
        self.c.path = './tests/data/out02391.dat'
        self.star_descriptions = [self.a, self.b, self.c]


    def test_read_lightcurve_sds(self):
        sds = reading.read_lightcurve_sds(self.star_descriptions)
        self.assertEqual(3, len(sds))


    @staticmethod
    def stardesc(id, ra, dec):
        return StarDescription(local_id=id,
                               coords=SkyCoord(ra, dec, unit='deg'))


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

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


    def test_extract_reference_frame_rotation(self):
        result = reading.extract_reference_frame_rotation('./tests/data/', 'WWCrA#30V_000000287_FLAT.fit')
        self.assertEqual(0.0, result)
        result = reading.extract_reference_frame_rotation('./tests/data/', 'WWCrA#30V_000000303_FLAT.fit')
        self.assertEqual(0.031, result)
        result = reading.extract_reference_frame_rotation('./tests/data/', 'WWCrA#30V_000000323_FLAT.fit')
        self.assertEqual(180.238, result)


    def test_extract_frame_from_summary_helper(self):
        result = reading.extract_frame_from_summary_helper('./tests/data/', "Ref.  image")
        self.assertEqual(('2458836.58742', '19.12.2019', '02:05:23',
                          '../data/ASASSN-V_J060000.76-310027.83/cleaned/2019/ASASSN-V_060000.76-310027.83#60V_000783664_FLAT.fit'), result)

    @staticmethod
    def stardesc(id, ra, dec):
        return StarDescription(local_id=id,
                               coords=SkyCoord(ra, dec, unit='deg'))


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

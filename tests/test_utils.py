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

logging.getLogger().setLevel(logging.DEBUG)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")

test_file_path = PurePath(os.getcwd(), 'tests', 'data')


class TestUtils(unittest.TestCase):

    def setUp(self) -> None:
        self.a = self.stardesc(1, 1, 1)
        self.b = self.stardesc(2, 2, 2)
        self.c = self.stardesc(3, 3, 3)
        self.star_descriptions = [self.a, self.b, self.c]
        self.set_catalog(self.star_descriptions[0], catalog_name="VSX")
        self.set_catalog(self.star_descriptions[0], catalog_name="CANDIDATES")
        self.set_catalog(self.star_descriptions[1], catalog_name="CANDIDATES")
        self.set_catalog(self.star_descriptions[2], catalog_name="VSX")


    def test_metadata_sorter(self):

        self.assertEqual(15, utils.metadata_sorter.get_string_number_part_or_default("oenutheo15"))
        self.assertEqual(1, utils.metadata_sorter.get_string_number_part_or_default("RMH-HMB-1"))
        self.assertEqual(7, utils.metadata_sorter.get_string_number_part_or_default("RMH-NEW-7"))
        self.assertEqual(sys.maxsize, utils.metadata_sorter.get_string_number_part_or_default("RMH-NEW"))

        stars = []
        for idx in range(1, 101):
            curr_star = self.stardesc(idx, random.random() * 360, random.random() * 90)
            curr_star.metadata = star_metadata.CatalogData(key="RMH-HMB", catalog_id=f'RMH-HMB-{idx}',
                                                           name=f'RMH-HMB-{idx}', coords=curr_star.coords, separation=0)
            stars.append(curr_star)
        random.shuffle(stars)
        result = utils.metadata_sorter(stars, metadata_id='RMH-HMB', name_variable='name')
        self.assertEqual("RMH-HMB-1", result[0].get_metadata("RMH-HMB").name)
        self.assertEqual("RMH-HMB-10", result[9].get_metadata("RMH-HMB").name)
        self.assertEqual("RMH-HMB-100", result[99].get_metadata("RMH-HMB").name)

        stars = []
        for idx in range(1, 101):
            curr_star = self.stardesc(idx, random.random() * 360, random.random() * 90)
            curr_star.metadata = star_metadata.CatalogData(key="RMH-HMB", catalog_id='', name='',
                                                           coords=curr_star.coords, separation=0)
            stars.append(curr_star)
        utils.metadata_sorter(stars, metadata_id="DOESNOTEXIST")  # expect no exception
        utils.metadata_sorter(stars, metadata_id="RMH-HMB")  # expect no exception

        stars = []
        for idx in range(1, 101):
            curr_star = self.stardesc(idx, random.random() * 360, random.random() * 90)
            curr_star.metadata = star_metadata.CatalogData(key="RMH-HMB", catalog_id=f'RMH-NEW{idx}',
                                                           name=f'RMH-NEW{idx}', coords=curr_star.coords, separation=0)
            stars.append(curr_star)
        random.shuffle(stars)
        result = utils.metadata_sorter(stars, metadata_id='RMH-HMB')
        self.assertEqual("RMH-NEW1", result[0].get_metadata("RMH-HMB").name)
        self.assertEqual("RMH-NEW10", result[9].get_metadata("RMH-HMB").name)
        self.assertEqual("RMH-NEW100", result[99].get_metadata("RMH-HMB").name)

        # MIXED stars
        stars = []
        random.seed(42)
        for idx in range(1, 101):
            curr_star = self.stardesc(idx, random.random() * 360, random.random() * 90)
            starname = f"RMH-HMB-{idx}" if random.random() < 0.5 else idx
            curr_star.metadata = star_metadata.\
                CatalogData(key="RMH-HMB", catalog_id=starname, name=starname,
                            coords=curr_star.coords, separation=0)
            stars.append(curr_star)
        stars[0].get_metadata("RMH-HMB").name="RMH-HMB"
        result = utils.metadata_sorter(stars, metadata_id="RMH-HMB")  # expect no exception
        self.assertEqual(2, result[0].get_metadata("RMH-HMB").name)
        self.assertEqual(91, result[47].get_metadata("RMH-HMB").name)
        self.assertEqual("RMH-HMB-3", result[49].get_metadata("RMH-HMB").name)
        self.assertEqual("RMH-HMB-15", result[53].get_metadata("RMH-HMB").name)


    def test_reject_outliers_iqr(self):
        np.random.seed(42)
        data = pd.DataFrame(np.random.random_sample(size=(100, 2)), columns=list('AB'))
        column = 'A'
        df = utils.reject_outliers_iqr(data, column, cut=30)
        self.assertEqual(len(data), len(df))
        df['A'][50] = -0.4
        df = utils.reject_outliers_iqr(df, column, cut=30)
        self.assertEqual(len(data) - 1, len(df))


    def test_concat_sd_lists(self):
        list1 = [self.a, self.b]
        list2 = [self.b, self.c]
        list3 = [self.c, self.b, self.a]
        result = utils.concat_sd_lists(list1, list2)
        self.assertEqual(3, len(result))
        result = utils.concat_sd_lists([], [])
        self.assertEqual(0, len(result))
        result = utils.concat_sd_lists(list1, list1)
        self.assertEqual(2, len(result))
        result = utils.concat_sd_lists(list1, list2, list3)
        self.assertEqual(3, len(result))


    def test_get_stars_with_metadata(self):
        vsx_descr = utils.get_stars_with_metadata(self.star_descriptions, catalog_name="CANDIDATES", exclude=["VSX"])
        self.assertEqual(1, len(vsx_descr))
        vsx_descr = utils.get_stars_with_metadata(self.star_descriptions, catalog_name="VSX", exclude=["CANDIDATES"])
        self.assertEqual(1, len(vsx_descr))
        vsx_descr = utils.get_stars_with_metadata(self.star_descriptions, catalog_name="VSX", exclude=["VSX",
                                                                                                       "CANDIDATES"])
        self.assertEqual(0, len(vsx_descr))


    def test_metadata_filter(self):
        # # Does this star have a catalog with catalog_name? Used in combination with filter()
        # def metadata_filter(star: StarDescription, catalog_name, exclude=[]):
        #     catalogs = star.get_metadata_list()
        #     return catalog_name in catalogs and len([x for x in exclude if x in catalogs]) == 0

        result = utils.metadata_filter(self.star_descriptions[0], catalog_name="CANDIDATES", exclude=["VSX"])
        self.assertFalse(result)
        result = utils.metadata_filter(self.star_descriptions[1], catalog_name="CANDIDATES", exclude=["VSX"])
        self.assertTrue(result)
        result = utils.metadata_filter(self.star_descriptions[1], catalog_name="CANDIDATES", exclude=[])
        self.assertTrue(result)
        result = utils.metadata_filter(self.star_descriptions[2], catalog_name="VSX", exclude=[])
        self.assertTrue(result)
        try:
            result = utils.metadata_filter(self.star_descriptions[2], catalog_name="VSX", exclude="VSX")
            self.fail("Should give an exception because exclude is not a list")
        except:
            pass


    @staticmethod
    def stardesc(id, ra, dec):
        return StarDescription(local_id=id,
                               coords=SkyCoord(ra, dec, unit='deg'))


    @staticmethod
    def set_catalog(star, catalog_name: str):
        star.metadata = CatalogData(key=catalog_name, catalog_id=f"{catalog_name}-{star.local_id}",
                                    name=f"{catalog_name}-{star.local_id}", separation=None, coords=star.coords,
                                    extradata=None)
        return star


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

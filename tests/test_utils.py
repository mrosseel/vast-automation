# from .context import src
import unittest
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

logging.getLogger().setLevel(logging.DEBUG)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")

test_file_path = PurePath(os.getcwd(), 'tests', 'data')


class TestUtils(unittest.TestCase):

    def test_sort_rmh_hmb(self):

        self.assertEqual(15, utils.metadata_sorter.get_metadata_name_number_part("oenutheo15"))
        self.assertEqual(1, utils.metadata_sorter.get_metadata_name_number_part("RMH-HMB-1"))
        self.assertEqual(7, utils.metadata_sorter.get_metadata_name_number_part("RMH-NEW-7"))

        stars = []
        for idx in range(1, 101):
            curr_star = self.stardesc(idx, random.random() * 360, random.random() * 90)
            curr_star.metadata = star_metadata.CatalogData(key="RMH-HMB", catalog_id=f'RMH-HMB-{idx}',
                                                           name=f'RMH-HMB-{idx}', coords=curr_star.coords, separation=0)
            stars.append(curr_star)
        random.shuffle(stars)
        name_extractor = lambda x: x.name
        result = utils.metadata_sorter(stars, metadata_id='RMH-HMB', name_extract=name_extractor)
        self.assertEqual("RMH-HMB-1", result[0].get_metadata("RMH-HMB").name)
        self.assertEqual("RMH-HMB-10", result[9].get_metadata("RMH-HMB").name)
        self.assertEqual("RMH-HMB-100", result[99].get_metadata("RMH-HMB").name)

        stars = []
        for idx in range(1, 101):
            curr_star = self.stardesc(idx, random.random() * 360, random.random() * 90)
            curr_star.metadata = star_metadata.CatalogData(key="RMH-HMB", catalog_id='', name='',
                                                           coords=curr_star.coords, separation=0)
            stars.append(curr_star)
        utils.metadata_sorter(stars)  # expect no exception

        stars = []
        for idx in range(1, 101):
            curr_star = self.stardesc(idx, random.random() * 360, random.random() * 90)
            curr_star.metadata = star_metadata.CatalogData(key="RMH-HMB", catalog_id=f'RMH-NEW{idx}',
                                                           name=f'RMH-NEW{idx}', coords=curr_star.coords, separation=0)
            stars.append(curr_star)
        random.shuffle(stars)
        result = utils.metadata_sorter(stars, metadata_id='RMH-HMB', name_extract=name_extractor)
        self.assertEqual("RMH-NEW1", result[0].get_metadata("RMH-HMB").name)
        self.assertEqual("RMH-NEW10", result[9].get_metadata("RMH-HMB").name)
        self.assertEqual("RMH-NEW100", result[99].get_metadata("RMH-HMB").name)


    def test_reject_outliers_iqr(self):
        np.random.seed(42)
        data = pd.DataFrame(np.random.random_sample(size=(100, 2)), columns=list('AB'))
        column = 'A'
        df = utils.reject_outliers_iqr(data, column, cut=30)
        self.assertEqual(len(data), len(df))
        df['A'][50] = -0.4
        df = utils.reject_outliers_iqr(df, column, cut=30)
        self.assertEqual(len(data) - 1, len(df))


    @staticmethod
    def stardesc(id, ra, dec):
        return StarDescription(local_id=id,
                               coords=SkyCoord(ra, dec, unit='deg'))


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

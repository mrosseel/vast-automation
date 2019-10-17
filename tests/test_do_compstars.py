# from .context import src
import unittest
import do_compstars
from comparison_stars import ComparisonStars
from star_description import StarDescription
from astropy.coordinates import SkyCoord
import logging
import main_vast
import os
from pathlib import PurePath
import logging
from pandas import DataFrame

logging.getLogger().setLevel(logging.DEBUG)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")

test_file_path = PurePath(os.getcwd(), 'tests', 'data')


class TestDoCompstars(unittest.TestCase):

    def test_calculate_mean_value_ensemble_photometry(self):
        data = {'JD':  ['1'],
                'Vrel': [15.414],
                'err': [0.012]
                }
        df = DataFrame(data, columns=['JD', 'Vrel', 'err'])
        # {JD, (mag, magerr)}
        observations_1 = {'1': (11.775, 0.001)}
        observations_2 = {'1': (12.220, 0.0012)}
        observations_3 = {'1': (13.114, 0.0021)}
        observations = [observations_1, observations_2, observations_3]
        catalogmags = [11.8, 12.2, 13.1]
        comp_stars = ComparisonStars(None, None, observations, catalogmags, None)

        realV, realErr = do_compstars.calculate_mean_value_ensemble_photometry(df, comp_stars)
        self.assertEqual(15.411, realV[0])
        self.assertEqual(0.0121, round(realErr[0], 4))


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

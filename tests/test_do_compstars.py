# from .context import src
import unittest
import do_compstars
from comparison_stars import ComparisonStars
from star_description import StarDescription
from astropy.coordinates import SkyCoord
import logging
import main_vast
import os
from pathlib import Path
import logging
from pandas import DataFrame
import utils

logging.getLogger().setLevel(logging.DEBUG)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")

test_file_path = Path(os.getcwd(), 'tests', 'data')


class TestDoCompstars(unittest.TestCase):

    def test_calculate_mean_value_ensemble_photometry(self):
        data = {'JD': ['1'],
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

        realV, realErr = do_compstars.calculate_ensemble_photometry(df, comp_stars,
                                                                    do_compstars.mean_value_ensemble_method)
        self.assertEqual(15.411, realV[0])
        self.assertEqual(0.0121, round(realErr[0], 4))

    def test_calculate_weighted_value_ensemble_photometry(self):
        data = {'JD': ['1'],
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

        realV, realErr = do_compstars.calculate_ensemble_photometry(df, comp_stars,
                                                                    do_compstars.weighted_value_ensemble_method)
        self.assertEqual(15.415, round(realV[0], 3))
        self.assertEqual(0.0121, round(realErr[0], 4))


    def test_get_calculated_compstars(self):
        stars = [self.stardesc(1, 1, 1, 10, 0.01, 10),
                 self.stardesc(4283, 3, 3, 12, 0.02, 10),
                 self.stardesc(132, 10.24496, 9.96736, 13, 0.02, 10),
                 self.stardesc(2, 10.24490, 9.96730, 12, 0.01, 10),
                 self.stardesc(3, 10.24490, 9.96730, 12, 0.01, 10)]  # not there
        stardict = utils.get_localid_to_sd_dict(stars)
        # def get_calculated_compstars(vastdir, stardict: StarDict, ref_jd, maglimit=15, starlimit=1000):
        ids, sds = do_compstars.get_calculated_compstars(test_file_path, stardict, ref_jd=1)
        self.assertEqual(4, len(ids))


    def test_closest_compstar_ids(self):
        stars = [self.stardesc(1, 1, 1, 10, 0.01, 10),
                 self.stardesc(4283, 3, 3, 12, 0.02, 10),
                 self.stardesc(132, 10.24496, 9.96736, 13, 0.02, 10),
                 self.stardesc(2, 10.24490, 9.96730, 12, 0.01, 10),
                 self.stardesc(3, 10.24490, 9.96730, 12, 0.01, 10)]  # not there
        result = do_compstars._closest_compstar_ids(stars[0],
                                                    ComparisonStars([x.local_id for x in stars[1:]],
                                                                   [x for x in stars[1:]], None, None, None))
        self.assertEqual([4283, 2, 3, 132], result)


    def test_get_list_of_likely_constant_stars(self):
        result = do_compstars._get_list_of_likely_constant_stars(test_file_path)
        self.assertEqual(6609, len(result))


    def stardesc(self, id, ra, dec, vmag, e_vmag, obs):
        return StarDescription(local_id=id,
                               coords=SkyCoord(ra, dec, unit='deg'), vmag=vmag, e_vmag=e_vmag, obs=obs)


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()


# from .context import src
import unittest
import do_calibration
import reading
import star_metadata
import utils
from star_description import StarDescription
from astropy.coordinates import SkyCoord
import logging
import main_vast
import os
from pathlib import PurePath
import logging
logging.getLogger().setLevel(logging.DEBUG)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")

test_file_path = PurePath(os.getcwd(), 'tests', 'data')


class TestMainVast(unittest.TestCase):

    def test_tag_candidates(self):
        stars = [self.stardesc(1, 1, 1), self.stardesc(2, 3, 3), self.stardesc(3129, 10.24496, 9.96736),
                 self.stardesc(5711, 10.24490, 9.96730)]
        main_vast.tag_candidates(test_file_path, stars)
        test = utils.get_stars_with_metadata(stars, "CANDIDATE")
        self.assertEqual(2, len(test))

    def test_tag_selected(self):
        # 3154, 515, 2269, 2635
        stars = [self.stardesc(515, 12.7725921, 12.0036631), self.stardesc(2269, 12.2686600, 12.0382637),
                 self.stardesc(3154, 12.1520503, 12.0935881),
                 self.stardesc(2635, 10.24490, 9.96730), self.stardesc(1, 10, 9)]
        main_vast.read_and_tag_localid(PurePath(test_file_path, "wwcra2015_starlist.txt"),
                                       utils.get_localid_to_sd_dict(stars))
        test = utils.get_stars_with_metadata(stars, "SITE")
        self.assertEqual(4, len(test))
        test = utils.get_stars_with_metadata(stars, "SELECTEDFILE")
        self.assertEqual(4, len(test))


    def stardesc(self, id, ra, dec):
        return StarDescription(local_id=id,
                               coords=SkyCoord(ra, dec, unit='deg'))


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

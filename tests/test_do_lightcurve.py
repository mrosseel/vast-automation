from .context import src
import unittest
import do_calibration
from star_description import StarDescription
from astropy.coordinates import SkyCoord
import logging

class TestDoCalibration(unittest.TestCase):

    def test_synthetic_c(self):
        stars = []
        stars.append(self.stardesc(1, 1, 1))
        stars.append(self.stardesc(1, 3, 3))
        stars.append(self.stardesc(1, 10.24496,  9.96736))
        stars.append(self.stardesc(1, 10.24490,  9.96730))
        result = do_calibration.add_vsx_names_to_star_descriptions(stars, 'FIXME_SETTINGS_NOT_WORKING', max_separation=0.1)
        for entry in result:
            print(entry, '\n')
        total = 0
        for entry in result:
            if entry.match:
                total += 1
        self.assertEqual(1, total)


    def stardesc(self, id, ra, dec):
        return StarDescription(local_id=id,
                               coords=SkyCoord(ra, dec, unit='deg'))


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

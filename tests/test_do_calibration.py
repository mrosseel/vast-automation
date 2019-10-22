import unittest
import do_calibration
from star_description import StarDescription
from astropy.coordinates import SkyCoord
import logging

class TestDoCalibration(unittest.TestCase):

    def test_add_info_to_star_description(self):
        stars = [self.stardesc(1, 1, 1), self.stardesc(1, 3, 3), self.stardesc(1, 10.24496, 9.96736),
                 self.stardesc(1, 10.24490, 9.96730)]
        star0_id = id(stars[0])
        result = do_calibration.add_info_to_star_description(stars[0], 10, 0.1, "BLA", "BLA",
                                                             SkyCoord(ra=8, dec=5, unit="deg"))
        self.assertEqual(star0_id, id(result))
        self.assertEqual(stars[0].vmag, 10)
        self.assertEqual(stars[0].e_vmag, 0.1)


    def test_add_VSX(self):
        stars = [self.stardesc(1, 1, 1), self.stardesc(1, 3, 3), self.stardesc(1, 10.24496, 9.96736),
                 self.stardesc(1, 10.24490, 9.96730)]
        result = do_calibration.add_vsx_names_to_star_descriptions(stars, 'FIXME', max_separation=0.1)
        for entry in result:
            print(entry, '\n')
        total = 0
        for entry in result:
            if entry.metadata:
                total += 1
        self.assertEqual(1, total)


    def stardesc(self, id, ra, dec):
        return StarDescription(local_id=id,
                               coords=SkyCoord(ra, dec, unit='deg'))


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

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
        # matching vsx star is the ninth from the vsx_mini.bin
        # >>> import vsx_pickle
        # >>> dicty = vsx_pickle.read('../tests/data/vsx_mini.bin')
        # >>> for entry in dicty:
        #     ...     print(dicty[entry][9])
        # 0.025
        # -59.74675
        # {'id': 9, 'OID': 170899, 'Name': 'UNSW-V 312', 'Type': 'EA', 'l_Period': nan,
        # 'Period': 1.05762, 'u_Period': nan}

        stars = [self.stardesc(1, 1, 1), self.stardesc(1, 3, 3), self.stardesc(1, 10.24496, 9.96736),
                 self.stardesc(1, 0.025, -59.74675)]
        sds, result_ids = do_calibration.add_vsx_names_to_star_descriptions(stars, './tests/data/vsx_mini.bin',
                                                                            max_separation=0.1)
        for entry in sds:
            print(entry, '\n')
        total = 0
        for entry in sds:
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

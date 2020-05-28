# from .context import src
import unittest
import do_compstars
from star_description import StarDescription
from ucac4 import UCAC4
from comparison_stars import ComparisonStars
import os
from pathlib import PurePath
import logging
from pandas import DataFrame
from pathlib import Path
from astropy.coordinates import SkyCoord

test_file_path = PurePath(os.getcwd(), 'tests', 'data')
epsilon = 0.0000001


class TestUcac4(unittest.TestCase):
    def setUp(self) -> None:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
        self.UCAC4_ID = 'UCAC4 001-000003'
        # self.ucac4 = UCAC4(ucac_path=Path('tests/data/'))
        self.ucac4 = UCAC4()


    def test_get_ucactuple_for_zone_and_runnrs(self):
        # ucac4.get_ucac4_star_description_raw('UCAC4 001-000003')
        result = self.ucac4.get_ucactuples_for_zone_and_runnrs(1, [3])
        self.assertEqual(result[0][0].ra, 18290451)


    def test_get_ucactuple_from_id(self):
        result = self.ucac4.get_ucactuple_from_id(self.UCAC4_ID)
        self.assertEqual(result[0].ra, 18290451)


    def test_get_ucac4_star_description(self):
        result = self.ucac4.get_star_description_from_id(self.UCAC4_ID)
        logging.debug(f"test_get_ucac4_star_description result: {result}")
        #     local_id: None, aavso_id: UCAC4 001-000003, coords: <SkyCoord (ICRS): (ra, dec) in deg
        # (5.08068083, -89.81054306)>, xy: None, None, vmag: 10.986, nr matches: 0, matches: [], path: None
        self.assertEqual(result.aavso_id, self.UCAC4_ID)
        self.assertTrue(result.coords.ra.deg - 5.08068083 < epsilon)
        self.assertTrue(result.coords.dec.deg - -89.81054306 < epsilon)
        self.assertEqual(result.vmag, 10.986)


    def test_name_to_zone_and_run_nr(self):
        zone, runnr = self.ucac4.ucac_id_to_zone_and_run_nr(self.UCAC4_ID)
        self.assertEqual(zone, 1)
        self.assertEqual(runnr, 3)


    def test_zone_and_run_nr_to_name(self):
        ucac4_id = self.ucac4.zone_and_run_nr_to_name(1, 3)
        self.assertEqual(ucac4_id, self.UCAC4_ID)


    def test_get_ucac4_id(self):
        # Compstar match: UCAC4 231-154752 with 3174 (271.2347344444444 deg, -43.84581611111111 deg)
        # Compstar match: UCAC4 232-147677 with 2620 (271.2807819444444 deg, -43.77729194444444 deg)
        ra, dec = 271.23473, -43.845816
        target = SkyCoord(ra, dec, unit='deg')
        result = self.ucac4.get_sd_from_ra_dec(ra, dec)
        self.assertEqual("UCAC4 231-154752", result.aavso_id)
        self.assertEqual(12.107, result.vmag)
        print("diff is ", target.separation(result.coords))
        ra, dec = 271.28078, -43.77729
        result = self.ucac4.get_sd_from_ra_dec(ra, dec)
        self.assertEqual("UCAC4 232-147677", result.aavso_id)
        self.assertEqual(12.314, result.vmag)
        print("diff is ", target.separation(result.coords))
        # WARNING Did not find a UCAC4 match for 274.26921036101504, -88.73827035367276, 0.02.
        # Buckets: range(1098, 1099), zones: [7],smallest dist: 1000
        # ra, dec = 274.26921036101504, -88.73827035367276
        ra, dec = 274.2390118, -88.7390245
        result = self.ucac4.get_sd_from_ra_dec(ra, dec, tolerance_deg=0.19)

        self.assertEqual("UCAC4 007-002358", result.aavso_id)
        # line 9738 (1 based)
        # 2358    0   7 1098
        # 2358    0   7 1099
        # 2358    2   7 1100
        self.assertEqual(20.0, result.vmag)
        print("diff is ", target.separation(result.coords))

        ra, dec = 271.0133100, -43.5485520
        # target = SkyCoord(ra, dec, unit='deg')
        result = self.ucac4.get_sd_from_ra_dec(ra, dec)
        self.assertEqual("UCAC4 233-155696", result.aavso_id)


    def test_get_zone_for_dec(self):
        self.assertEqual(1, self.ucac4.get_zone_for_dec(-90))
        self.assertEqual(2, self.ucac4.get_zone_for_dec(-89.8))
        self.assertEqual(3, self.ucac4.get_zone_for_dec(-89.6))
        self.assertEqual(4, self.ucac4.get_zone_for_dec(-89.4))
        self.assertEqual(5, self.ucac4.get_zone_for_dec(-89.2))
        self.assertEqual(900, self.ucac4.get_zone_for_dec(90))


    def test_get_ra_bucket(self):
        self.assertEqual(1440, self.ucac4.ra_bin_index(360))  # 1440 but zero based
        self.assertEqual(1, self.ucac4.ra_bin_index(0))
        self.assertEqual(720, self.ucac4.ra_bin_index(180))  # 720 but zero based


    def test_mag_error(self):
        # 'UCAC4 232-146904'
        # <SkyCoord (ICRS): (ra, dec) in deg (270.91140731, -43.64018836)>
        # sd: StarDescription = self.ucac4.get_ucac4_sd(270.91140731, -43.64018836)
        raw = self.ucac4.get_ucactuples_for_zone_and_runnrs(232, [146904])
        # Vmag	11.957 	mag, e_Vmag	01 cmag
        self.assertEqual(-1, raw[0][0].apass_mag_sigma_V)
        sd: StarDescription = self.ucac4.get_star_description_from_tuple(raw[0])
        self.assertEqual(.01, sd.vmag_err)

        # first entry of out.sam
        # raw = self.ucac4.get_ucac4_details_raw('451', [133336])
        # self.assertEqual(5, raw[0][0].apass_mag_sigma_V)
        # sd: StarDescription = self.ucac4.get_ucac4_star_description_fromtuple(*raw[0])
        # self.assertEqual(.05, sd.e_vmag)
        # print("bl")


    def test_get_ucac4_range_tuples(self):
        result = self.ucac4.get_region_minimal_star_tuples(40, 40, 1)
        self.assertEqual(2480, len(result))


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

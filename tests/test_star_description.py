import unittest
from star_description import StarDescription
from astropy.coordinates import SkyCoord
import logging
from star_metadata import CompStarData

class TestStarDescription(unittest.TestCase):

    def test_metadata(self):
        star = StarDescription(local_id=id, coords=SkyCoord(10, 11, unit='deg'))
        self.assertEqual([], star.metadata)
        star.metadata = CompStarData()
        self.assertEqual(1, len(star.metadata))
        star.metadata = CompStarData()
        self.assertEqual(2, len(star.metadata))



if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

import unittest
from star_description import StarDescription
from astropy.coordinates import SkyCoord
import logging
from star_metadata import CompStarData
from star_metadata import StarMetaData


class TestStarDescription(unittest.TestCase):

    def test_metadata(self):
        star = StarDescription(local_id=id, coords=SkyCoord(10, 11, unit='deg'))
        self.assertEqual({}, star.metadata)
        star.metadata = CompStarData([1])
        self.assertEqual(1, len(star.metadata))
        star.metadata = StarMetaData()
        self.assertEqual(2, len(star.metadata))
        try:
            star.set_metadata(StarMetaData(), True)
        except:
            return
        self.fail("Adding duplicate keys to the dictionary should fail")



if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

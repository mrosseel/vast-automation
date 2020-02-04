import unittest
import do_charts_vast
import numpy as np
from star_description import StarDescription
from astropy.coordinates import SkyCoord
import logging


class TestDoChartsVast(unittest.TestCase):

    def test_shift_to_epoch(self):
        np.random.seed(42)
        t_np = np.random.uniform(2456810, 2456820, 10)
        t_str = np.array(list(map(str, t_np)))
        epoch = '2456810.5808361215'
        t_np_zeroed = do_charts_vast.shift_to_epoch(epoch, t_np, t_str)
        self.assertEqual(10, len(t_np_zeroed))
        self.assertEqual(0, t_np_zeroed[6])


    def stardesc(self, id, ra, dec):
        return StarDescription(local_id=id,
                               coords=SkyCoord(ra, dec, unit='deg'))


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

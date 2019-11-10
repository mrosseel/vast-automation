# from .context import src
import unittest
import do_calibration
import star_metadata
import utils
from star_description import StarDescription
from astropy.coordinates import SkyCoord
import logging
import main_vast
import os
from pathlib import PurePath
import logging
import reading
import do_aavso_report
import aavso
import pandas as pd
from timeit import default_timer as timer

logging.getLogger().setLevel(logging.DEBUG)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")

test_file_path = PurePath(os.getcwd(), 'tests', 'data')


class TestDoAAVSO(unittest.TestCase):

    def stardesc(self, id, ra, dec):
        return StarDescription(local_id=id,
                               coords=SkyCoord(ra, dec, unit='deg'))


    def test_do_aavso(self):
        df = reading.read_lightcurve_vast('./tests/data/out02391.dat')
        print(df.head())
        sd = self.stardesc(2391, 10, 10)
        # do_aavso_report.report(sd, df)
        data = {
            'name': 'var_display_name',
            'date': '212312321.12',
            'magnitude': 12.0,
            'magnitude_error': 0.01,
            'filter': 'V',
            'transformed': 'YES',
            'magnitude_type': 'STD',
            'comparison_name': 'ENSEMBLE',
            'comparison_magnitude': 'na',
            'check_name': 'check_display_name',
            'check_magnitude': 'comparison_star_vmag',
            'airmass': 'airmass',
            'group': 'na',
            'chart': 'na',
            'notes': 'na'
        }
        start = timer()
        with open('./tests/data/aavso_out.txt', 'w') as fp:
            writer = aavso.ExtendedFormatWriter(fp, 'RMH', software='munipack-automation', type='EXTENDED',
                                                obstype='CCD')
            for i in range(100000):
                writer.addrow(data)
            s = "".join(writer.data)
            writer.flush()
        end = timer()
        import hashlib
        m = hashlib.sha256()
        m.update(s.encode('utf-8'))
        result = m.hexdigest()
        self.assertEqual("abe0456d54f6070739fec226f25f24b4790d756f1cd7c212d0bfcd3eb8dc52ff", result)
        print("timing is", end-start)
    # print(df.describe())
        # print(df.info())


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

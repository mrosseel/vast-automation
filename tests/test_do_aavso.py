# from .context import src
import unittest
import do_calibration
import star_metadata
import utils
from comparison_stars import ComparisonStars
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
import numpy as np
from pandas import DataFrame
import pandas as pd
from timeit import default_timer as timer

logging.getLogger().setLevel(logging.DEBUG)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")

test_file_path = PurePath(os.getcwd(), 'tests', 'data')


def massage_df_for_phase_plot(df: DataFrame):
    df['floatJD'] = df['JD'].astype(np.float)
    df['realV'] = df['Vrel']
    df['realErr'] = df['err']
    return df


def stardesc(id, ra, dec):
    return StarDescription(local_id=id,
                           coords=SkyCoord(ra, dec, unit='deg'))


class TestDoAAVSO(unittest.TestCase):

    def setUp(self) -> None:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
        self.df = reading.read_lightcurve_vast('./tests/data/out02391.dat')
        self.sd = stardesc(2391, 10, 10)
        do_calibration.add_info_to_star_description(self.sd, 10, 0.1, "UCAC4-Numero2391", "UCAC4", self.sd.coords)
        self.sd_comp1 = stardesc(1, 10, 10)
        do_calibration.add_info_to_star_description(self.sd_comp1, 12.1, 0.1, "UCAC4-Numero1", "UCAC4", self.sd.coords)
        self.sd_comp5 = stardesc(5, 10, 10)
        do_calibration.add_info_to_star_description(self.sd_comp5, 11.1, 0.1, "UCAC4-Numero5", "UCAC4", self.sd.coords)
        self.sd.set_metadata(star_metadata.CompStarData(compstar_ids=[1, 5]))
        # ComparisonStars(comparison_stars_ids, comparison_stars_1_sds, comp_observations, comp_catalogmags,
        #                 comp_catalogerr)
        observations = [{'2457236.66302': 12.2}, {'2457236.66302': 11.2}]
        self.comp_stars = ComparisonStars([1, 5], [self.sd_comp1, self.sd_comp5], observations, [12, 11], [0.1, 0.1])


    def test_aavso(self):
        data = {
            'name': 'var_display_name',
            'date': '212312321.12',
            'magnitude': 12.0,
            'magnitude_error': 0.01,
            'filter': 'V',
            'transformed': 'NO',
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
        self.assertEqual("f042e25b9e6078d03b1e76e8ce9a12ac64ec1d16edbdff0fd452ca4a44ab2555", result)
        print("timing is", end - start)


    def test_do_aavso_report(self):
        # def report(star: StarDescription, df_curve: DataFrame, target_dir: Path, vastdir: str, sitelat, sitelong,
        #            sitealt, filter=None, observer='RMH', chunk_size=None):
        import toml
        settings = toml.load('./tests/data/testsettings.txt')
        do_aavso_report.report(self.sd, massage_df_for_phase_plot(self.df), self.comp_stars, './tests/data',
                               sitelat=settings['sitelat'], sitelong=settings['sitelong'], sitealt=settings['sitealt'],
                               observer=settings['observer'], camera_filter='V')


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

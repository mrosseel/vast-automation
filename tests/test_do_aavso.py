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

test_file_path = PurePath(os.getcwd(), "tests", "data")


def massage_df_for_phase_plot(df: DataFrame):
    df["floatJD"] = df["JD"].astype(np.float)
    df["realV"] = df["Vrel"]
    df["realErr"] = df["err"]
    return df


def stardesc(id, ra, dec):
    return StarDescription(local_id=id, coords=SkyCoord(ra, dec, unit="deg"))


class TestDoAAVSO(unittest.TestCase):
    def setUp(self) -> None:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
        self.df = reading.read_lightcurve_vast("./tests/data/out02391.dat")
        self.sd = stardesc(2391, 10, 10)
        do_calibration.add_catalog_data_to_sd(
            self.sd, 10, 0.1, "UCAC4-Numero2391", "UCAC4", self.sd.coords
        )
        self.sd_comp1 = stardesc(1, 10, 10)
        do_calibration.add_catalog_data_to_sd(
            self.sd_comp1, 12.1, 0.1, "UCAC4-Numero1", "UCAC4", self.sd.coords
        )
        self.sd_comp5 = stardesc(5, 10, 10)
        do_calibration.add_catalog_data_to_sd(
            self.sd_comp5, 11.1, 0.1, "UCAC4-Numero5", "UCAC4", self.sd.coords
        )
        self.sd.set_metadata(star_metadata.CompStarData(compstar_ids=[1, 5]))
        # ComparisonStars(comparison_stars_ids, comparison_stars_1_sds, comp_observations, comp_catalogmags,
        #                 comp_catalogerr)
        observations = [
            {"2457236.66302": (12.2, 0.01)},
            {"2457236.66302": (11.2, 0.01)},
        ]
        self.comp_stars = ComparisonStars(
            [1, 5], [self.sd_comp1, self.sd_comp5], observations, [12, 11], [0.1, 0.1]
        )
        import toml

        self.settings = toml.load("./tests/data/testsettings.txt")

    def test_aavso(self):
        data = {
            "name": "var_display_name",
            "date": "212312321.12",
            "magnitude": 12.0,
            "magnitude_error": 0.01,
            "filter": "V",
            "transformed": "NO",
            "magnitude_type": "STD",
            "comparison_name": "ENSEMBLE",
            "comparison_magnitude": "na",
            "check_name": "check_display_name",
            "check_magnitude": "comparison_star_vmag",
            "airmass": 2.5,
            "group": "na",
            "chart": "na",
            "notes": "na",
        }
        start = timer()
        try:
            aavso_out_txt = "./tests/data/aavso_out.txt"
            with open(aavso_out_txt, "w") as fp:
                writer = aavso.ExtendedFormatWriter(
                    fp,
                    "RMH",
                    software="munipack-automation",
                    type="EXTENDED",
                    obstype="CCD",
                    location=(
                        self.settings["sitelat"],
                        self.settings["sitelong"],
                        self.settings["sitealt"],
                    ),
                )
                for i in range(100):
                    writer.addrow(data)
                s = "".join(writer.data)
                writer.flush()
        except:
            print("exception")
        end = timer()
        # NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
        # var_display_name,212312321.12,12.000,0.010,V,NO,STD,ENSEMBLE,na,check_display_name,comparison_star_vmag,2.50,na,na,na
        df = reading.read_aavso_lightcurve(aavso_out_txt)
        self.assertEqual(df.iloc[0]["NAME"], data["name"])
        self.assertEqual(df.iloc[0]["DATE"], str(data["date"]))
        self.assertEqual(df.iloc[0]["MAG"], data["magnitude"])
        self.assertEqual(df.iloc[0]["MERR"], data["magnitude_error"])
        self.assertEqual(df.iloc[0]["FILT"], data["filter"])
        self.assertEqual(df.iloc[0]["TRANS"], data["transformed"])
        self.assertEqual(df.iloc[0]["MTYPE"], data["magnitude_type"])
        self.assertEqual(df.iloc[0]["CNAME"], data["comparison_name"])
        self.assertEqual(df.iloc[0]["CMAG"], data["comparison_magnitude"])
        self.assertEqual(df.iloc[0]["KNAME"], data["check_name"])
        self.assertEqual(df.iloc[0]["KMAG"], data["check_magnitude"])
        self.assertEqual(df.iloc[0]["AMASS"], data["airmass"])
        self.assertEqual(df.iloc[0]["GROUP"], data["group"])
        self.assertEqual(df.iloc[0]["CHART"], data["chart"])
        self.assertEqual(df.iloc[0]["NOTES"], data["notes"])

        print("timing is", end - start)

    def test_do_aavso_report(self):
        # def report(star: StarDescription, df_curve: DataFrame, target_dir: Path, vastdir: str, sitelat, sitelong,
        #            sitealt, filter=None, observer='RMH', chunk_size=None):
        do_aavso_report.report(
            self.sd,
            massage_df_for_phase_plot(self.df),
            self.comp_stars,
            "./tests/data",
            target_dir="./tests/data",
            sitelat=self.settings["sitelat"],
            sitelong=self.settings["sitelong"],
            sitealt=self.settings["sitealt"],
            observer=self.settings["observer"],
            camera_filter="V",
        )


if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    unittest.main()

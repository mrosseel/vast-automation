from typing import List, Dict, Tuple
import os
import pandas as pd
import numpy as np
import errno
import re
import glob
import logging
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from collections import namedtuple

import do_calibration
import utils
from utils import StarDict
from star_description import StarDescription


# - 1st column - JD(TT) (default) or JD(UTC) (if VaST was started with "-u" flag)
# - 2nd column - magnitude (with respect to the background level on the reference image if an
#                absolute calibration was not done yet)
# - 3rd column - estimated magnitude error
# - 4th column - X position of the star on the current frame (in pixels)
# - 5th column - Y position of the star on the current frame (in pixels)
# - 6th column - diameter of the circular aperture used to measure the current frame (in pixels)
# - 7th column - file path corresponding to the current frame
def read_lightcurve_vast(starpath: str):
    logging.debug(f"Read lightcurve at path {starpath}")
    return pd.read_csv(
        starpath,
        delim_whitespace=True,
        names=["JD", "Vrel", "err", "X", "Y", "aperture?", "file"],
        usecols=["JD", "Vrel", "err", "X", "Y", "aperture?", "file"],
        dtype={"JD": str},
    )


def read_lightcurve_ids(star_ids: List[int], stardict: StarDict):
    result = []
    for star_id in star_ids:
        df = read_lightcurve_vast(stardict[star_id].path)
        result.append(df)
    return result


def read_lightcurve_sds(sds: List[StarDescription]):
    return list(map(lambda x: read_lightcurve_vast(x.path), sds))


def read_aavso_lightcurve(aavso_file: str):
    return pd.read_csv(
        aavso_file,
        sep=",",
        header=None,
        index_col=False,
        comment="#",
        names=[
            "NAME",
            "DATE",
            "MAG",
            "MERR",
            "FILT",
            "TRANS",
            "MTYPE",
            "CNAME",
            "CMAG",
            "KNAME",
            "KMAG",
            "AMASS",
            "GROUP",
            "CHART",
            "NOTES",
        ],
        dtype={"DATE": str},
    )


def trash_and_recreate_dir(adir: str):
    os.system('rm -fr "%s"' % adir)
    # shutil.rmtree(dir, ignore_errors=True)
    create_dir(adir)


def create_dir(adir: str):
    os.makedirs(adir, exist_ok=True)


def reduce_star_list(star_list_1, the_path):
    """
    takes a star_list and a dir, and returns a reduced star list - all stars which already have a file in that dir are
    removed
    """
    the_dir = os.listdir(the_path)
    the_dir.sort()
    found = []
    for filename in the_dir:
        found.append(filename_to_star(filename))
    logging.info(f"Found {len(found)} stars already processed in {the_path}")
    return [item for item in star_list_1 if item not in found]


# takes a filename and extracts the star number from it
def filename_to_star(filename):
    m = re.search(r"\d+", filename)
    return int(m.group(0).lstrip("0"))


# read the world positions and return them in a dictionary
# returns {'name': [ra.deg, dec.deg ]}
def read_world_positions(the_path):
    the_dir = os.listdir(the_path)
    the_dir.sort()
    results = {}
    for (
        name
    ) in the_dir:  # 'file' is a builtin type, 'name' is a less-ambiguous variable name.
        try:
            with open(
                the_path + name
            ) as f:  # No need to specify 'r': this is the default.
                results[filename_to_star(name)] = f.readlines()[0].split(" ")
        except IOError as exc:
            if (
                exc.errno != errno.EISDIR
            ):  # Do not fail if a directory is found, just ignore it.
                raise  # Propagate other kinds of IOError.
    return results


# get all possible stars with dict: JD, (mag, error)
def read_magdict_for_star(vastdir, star_id):
    stardict = {}
    starfile = Path(vastdir, star_to_dat(star_id))
    with open(starfile) as file:
        for line in file:
            splitline = line.split()
            # {JD, (mag, magerr)}
            stardict[str(splitline[0])] = (float(splitline[1]), float(splitline[2]))
    return stardict


#  blank_data is false if no background needs to be plotted, in that case all zeros are used as data
def get_fits_data(
    fits_file: str, blank_data: bool = False
) -> Tuple[List[float], float, float]:
    hdulist = fits.open(fits_file)
    data = hdulist[0].data.astype(float)
    shapex, shapey = hdulist[0].shape
    if blank_data:
        data = np.zeros(data.shape)
    return data, shapex, shapey


def read_vast_image_details_log(vastdir) -> pd.DataFrame:
    filename = Path(vastdir, "vast_image_details.log")
    with open(filename) as f:
        content = f.readline()
    max_len = len(content)
    col_specification = [
        (47, 60),
        (66, 70),
        (82, 89),
        (102, 107),
        (119, 124),
        (134, 138),
        (141, max_len),
    ]
    data = pd.read_fwf(
        filename,
        colspecs=col_specification,
        skiprows=0,
        names=("jd", "ap", "rotation", "detected", "matched", "status", "filename"),
        converters={
            "jd": float,
            "ap": float,
            "rotation": float,
            "detected": int,
            "matched": int,
            "status": str,
            "filename": str,
        },
    )
    return data


# TODO use the more modern 'read_vast_image_details_log' and filter out the ref frame
# get ref frame rotation from 'vast_image_details.log'
def extract_reference_frame_rotation(vastdir, reference_frame) -> float:
    filename = Path(vastdir, "vast_image_details.log")
    the_regex = re.compile(r"^.*rotation=\s*([0-9,.,-]+).*\s+(.+)$")
    with open(filename, "r") as infile:
        for line in infile:
            thesearch = the_regex.search(line)
            if thesearch and reference_frame in thesearch.group(2):
                return float(thesearch.group(1).strip())
    return 0.0


# get the mapping 'fits filename' -> rotation
def fitsfile_to_rotation_dict(vastdir) -> Dict[str, float]:
    df = read_vast_image_details_log(vastdir)
    rotation_dict = {}
    for index, row in df.iterrows():
        filename = Path(row.filename)
        rotation_dict[filename.name] = float(row.rotation)
    return rotation_dict


def jd_to_fitsfile_dict(vastdir) -> Dict[float, str]:
    df = read_vast_image_details_log(vastdir)
    jd_to_file_dict = {}
    for index, row in df.iterrows():
        filename = Path(row.filename)
        jd_to_file_dict[row.jd] = filename.name
    return jd_to_file_dict


# get the first image used by vast
def extract_first_frame(from_dir):
    return extract_frame_from_summary_helper(from_dir, "First image")


# get the reference image used by vast
def extract_reference_frame(from_dir):
    return extract_frame_from_summary_helper(from_dir, "Ref.  image")


# Get info from 'vast_summary.log'
# ref_jd, date, time, reference_frame
def extract_frame_from_summary_helper(from_dir, marker) -> Tuple[str, str, str, str]:
    # Ref.  image: 2458586.50154 13.04.2019 00:00:41   ../../inputfiles/TXCar/fits/TXCar#45V_000601040_FLAT.fit
    with open(Path(from_dir, "vast_summary.log")) as file:
        result = [
            re.findall(
                marker + r":\s+(\d+\.\d+)\s+([\d|\.]+)\s+([\d|\:]+)\s+(.*)", line
            )
            for line in file
        ]
    return [x for x in result if x != []][0][0]


# make a dict with mapping fits file -> image00007.cat
def extract_image_catalog(vastdir) -> float:
    filename = Path(vastdir, "vast_images_catalogs.log")
    the_regex = re.compile(r"^(.*) (.*)$")
    catalog_dict = {}
    with open(filename, "r") as infile:
        for line in infile:
            thesearch = the_regex.search(line)
            if thesearch:
                path = Path(thesearch.group(2))
                catalog_dict[path.name] = thesearch.group(1)
    return catalog_dict


# gets the nr of images used for photometry
def extract_images_used(from_dir):
    with open(from_dir + "vast_summary.log") as file:
        result = [re.findall(r"Images used for photometry (.*)", line) for line in file]
    return [x for x in result if x != []][0][0]


# Note: this file seems to give incorrect xy positions wrt reference frame
# Note: not used
# get all possible stars with their x/y position from a log file
# 14.460155 0.031190   215.230    19.626 out00007.dat
def read_data_m_sigma(vastdir) -> Dict[int, Tuple[int, int]]:
    stardict = {}
    PixelPos = namedtuple("PixelPos", "x y afile")
    with open(vastdir + "data.m_sigma") as file:
        for line in file:
            splitline = line.split()
            star_id = utils.get_starid_from_outfile(splitline[4])
            stardict[star_id] = PixelPos(
                float(splitline[2]), float(splitline[3]), splitline[4]
            )
    return stardict


ImageRecord = namedtuple("ImageRecord", "jd, x, y, file, rotation")


# for a certain star id, read the lightcurve and return for each observation the JD, X pos, Y pos and rotation
def get_star_jd_xy_rot(
    starid: int, vastdir: str
) -> Tuple[List[ImageRecord], Dict[str, float]]:
    lightcurvefile = f"out{starid:05}.dat"
    logging.info(f"file is {lightcurvefile}")
    df = read_lightcurve_vast(Path(vastdir, lightcurvefile))
    image_records = []
    rotation_dict = fitsfile_to_rotation_dict(vastdir)
    logging.info(f"rotation dict has {len(rotation_dict)} entries")
    for index, row in df.iterrows():
        filename = Path(row["file"]).name
        image_records.append(
            ImageRecord(
                float(row["JD"]),
                round(row["X"]),
                round(row["Y"]),
                filename,
                float(rotation_dict[filename]),
            )
        )
    return image_records, rotation_dict


def star_to_dat(star: int):
    return f"out{star:05}.dat"


def read_wcs_file(vastdir: str) -> Tuple[str, WCS]:
    wcs_file = Path(vastdir, "new-image.fits")
    # get wcs model from the reference header. Used in writing world positions and field charts
    try:
        wcs = do_calibration.get_wcs(wcs_file)
    except FileNotFoundError:
        wcs = None
    return wcs_file, wcs


# star id -> xpos, ypos, filename
StarPosDict = Dict[str, Tuple[float, float, str]]


# get a dict with star_id -> xpos ypos filename
def starid_to_xy_file_dict(vastdir: str) -> StarPosDict:
    stardict = {}
    PixelPos = namedtuple("PixelPos", "x y afile")
    with open(vastdir + "vast_list_of_all_stars.log") as file:
        for line in file:
            splitline = line.split()
            stardict[int(splitline[0])] = PixelPos(
                splitline[1], splitline[2], f"out{splitline[0]}.dat"
            )
    return stardict


# constructs a list of star descriptions with catalog matches according to args
def count_number_of_observations(vastdir):
    logging.info("Counting number of observations per star ...")
    obsdict = {}
    columns = [
        "Median magnitude",
        "idx00_STD",
        "X position of the star on the reference image [pix]",
        "Y position of the star on the reference image [pix]",
        "lightcurve file name",
        "idx01_wSTD",
        "idx02_skew",
        "idx03_kurt",
        "idx04_I",
        "idx05_J",
        "idx06_K",
        "idx07_L",
        "idx08_Npts",
        "idx09_MAD",
        "idx10_lag1",
        "idx11_RoMS",
        "idx12_rCh2",
        "idx13_Isgn",
        "idx14_Vp2p",
        "idx15_Jclp",
        "idx16_Lclp",
        "idx17_Jtim",
        "idx18_Ltim",
        "idx19_N3",
        "idx20_excr",
        "idx21_eta",
        "idx22_E_A",
        "idx23_S_B",
        "idx24_NXS",
        "idx25_IQR",
        "idx26_A01",
        "idx27_A02",
        "idx28_A03",
        "idx29_A04",
        "idx30_A05",
    ]
    df = pd.read_csv(
        Path(vastdir, "vast_lightcurve_statistics.log"),
        names=columns,
        delim_whitespace=True,
    )
    for index, row in df.iterrows():
        obsdict[row["lightcurve file name"]] = row["idx08_Npts"]
    return obsdict

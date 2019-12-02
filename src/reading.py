from typing import List
import os
import pandas as pd
import numpy as np
import errno
import re
import glob
import logging
import utils


# - 1st column - JD(TT) (default) or JD(UTC) (if VaST was started with "-u" flag)
# - 2nd column - magnitude (with respect to the background level on the reference image if an absolute calibration was not done yet)
# - 3rd column - estimated magnitude error
# - 4th column - X position of the star on the current frame (in pixels)
# - 5th column - Y position of the star on the current frame (in pixels)
# - 6th column - diameter of the circular aperture used to measure the current frame (in pixels)
# - 7th column - file path corresponding to the current frame
def read_lightcurve_vast(starpath: str):
    logging.debug(f"Read lightcurve at path {starpath}")
    return pd.read_csv(starpath, delim_whitespace=True,
                       names=['JD', 'Vrel', 'err', 'X', 'Y', 'unknown', 'file'],
                       usecols=['JD', 'Vrel', 'err', 'X', 'Y', 'unknown', 'file'], dtype={'JD': str})


def read_aavso_lightcurve(aavso_file: str):
    return pd.read_csv(aavso_file, sep=',', skiprows=7, header=None, index_col=False,
                       names=['NAME', 'DATE', 'MAG', 'MERR', 'FILT', 'TRANS', 'MTYPE', 'CNAME', 'CMAG', 'KNAME',
                              'KMAG', 'AMASS', 'GROUP', 'CHART,NOTES'])


def trash_and_recreate_dir(dir):
    os.system('rm -fr "%s"' % dir)
    # shutil.rmtree(dir, ignore_errors=True)
    create_dir(dir)


def create_dir(dir):
    os.makedirs(dir, exist_ok=True)


# takes a star_list and a dir, and returns a reduced star list - all stars which already have a file in that dir are removed
def reduce_star_list(star_list_1, the_path):
    the_dir = os.listdir(the_path)
    the_dir.sort()
    found = []
    for filename in the_dir:
        found.append(filename_to_star(filename))
    logging.info(f"Found {len(found)} stars already processed in {the_path}")
    return [item for item in star_list_1 if item not in found]


# takes a filename and extracts the star number from it
def filename_to_star(filename):
    m = re.search(r'\d+', filename)
    return int(m.group(0).lstrip('0'))


# read the world positions and return them in a dictionary
# returns {'name': [ra.deg, dec.deg ]}
def read_world_positions(the_path):
    the_dir = os.listdir(the_path)
    the_dir.sort()
    results = {}
    for name in the_dir:  # 'file' is a builtin type, 'name' is a less-ambiguous variable name.
        try:
            with open(the_path + name) as f:  # No need to specify 'r': this is the default.
                results[filename_to_star(name)] = f.readlines()[0].split(' ')
        except IOError as exc:
            if exc.errno != errno.EISDIR:  # Do not fail if a directory is found, just ignore it.
                raise  # Propagate other kinds of IOError.
    return results


# Select files conforming to the match_pattern using percentage which is between 0 and 1
def file_selector(the_dir, match_pattern, percentage=1) -> List[str]:
    matched_files = glob.glob(the_dir + match_pattern)
    desired_length = max(1, int(len(matched_files) * float(percentage)))
    logging.debug(f"Reading.file_selector: {the_dir + match_pattern}, "
                  f"total:{len(matched_files)}, desired:{desired_length}")
    np.random.seed(42)  # for the same percentage, we always get the same selection
    selected_files = np.random.choice(matched_files, size=desired_length, replace=False).tolist()
    return selected_files

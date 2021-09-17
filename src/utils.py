import glob
from functools import partial
from os import listdir
from os.path import isfile, join
import sys

import numpy as np
import star_description
from star_description import StarDescription, StarMetaData
from typing import List, Dict, Tuple
import multiprocessing as mp
from multiprocessing import cpu_count
import re
import logging
from astropy.coordinates import SkyCoord
from star_metadata import CatalogData
from collections import namedtuple
from pandas import DataFrame

# all info needed for ui purposes
StarUI = namedtuple(
    "StarInfo",
    "catalog_name, separation, extradata, filename_orig_no_ext, filename_no_ext, "
    "filename_orig_no_suff_no_ext, filename_no_suff_no_ext",
)
StarDict = Dict[int, StarDescription]


# Select files conforming to the match_pattern using percentage which is between 0 and 1
def file_selector(the_dir, match_pattern, percentage=1) -> List[str]:
    matched_files = glob.glob(the_dir + match_pattern)
    desired_length = max(1, int(len(matched_files) * float(percentage)))
    logging.debug(
        f"Reading.file_selector: {the_dir + match_pattern}, "
        f"total:{len(matched_files)}, desired:{desired_length}"
    )
    np.random.seed(42)  # for the same percentage, we always get the same selection
    selected_files = np.random.choice(
        matched_files, size=desired_length, replace=False
    ).tolist()
    return selected_files


def add_trailing_slash(the_path):
    return join(the_path, "")


def get_files_in_dir(mypath):
    return [f for f in listdir(mypath) if isfile(join(mypath, f))]


# out012345.dat -> 12345
def get_starid_from_outfile(outfile) -> int:
    m = re.search("out(.*).dat", outfile)
    try:
        result = int(m.group(1).lstrip("0"))
    except:
        logging.error(f"Could not extract the starid from file: {outfile}")
    return result


# returns a dict with the local_id as key
def get_localid_to_sd_dict(stars: List[StarDescription]) -> Dict[int, StarDescription]:
    cachedict = {}
    for sd in stars:
        cachedict[sd.local_id] = sd
    return cachedict


# filter a list of star descriptions on the presence of a catalog
def catalog_filter(star: StarDescription, catalog_name):
    return star.has_metadata(catalog_name)


# filters a DataFrame with a floatJD column according to julian dates
def jd_filter_df(df: DataFrame, jdfilter: List[float]):
    """ takes a list of 2 julian dates and uses these so the region between them is not used. The DataFrame needs a column named 'floatJD' """
    if jdfilter is not None:
        logging.debug(
            f"Before jd_filter_df(): jdfilter is: {jdfilter}, len is {len(df)}"
        )
        df = df[(df.floatJD <= jdfilter[0]) | (df.floatJD >= jdfilter[1])]
        logging.debug(f"After jd_filter_df(): len is {len(df)}")
        if len(df) < 2:
            logging.warning(
                f"Applying the jdfilter caused the lightcurve to contain less than 2 points! "
                f"Everything between {jdfilter[0]} and {jdfilter[1]} is thrown away"
            )
    return df


def jd_filter_array(jds, values, jdfilter: List[float]):
    """ takes a JD array and a value array (mags) together with a min/max jdfilter list """
    if jdfilter is not None:
        logging.debug(
            f"jd_filter_array(): jdfilter is: {jdfilter}, of type {type(jdfilter)}"
        )
        return zip(
            *filter(
                lambda jdsvaluezip: jdsvaluezip[0] <= jdfilter[0]
                or jdsvaluezip[0] >= jdfilter[1],
                zip(jds, values),
            )
        )
    else:
        return jds, values


def get_hms_dms(coord: SkyCoord):
    return "{:2.0f}h {:02.0f}m {:02.2f}s  {:2.0f}d {:02.0f}' {:02.2f}\"".format(
        coord.ra.hms.h,
        abs(coord.ra.hms.m),
        abs(coord.ra.hms.s),
        coord.dec.dms.d,
        abs(coord.dec.dms.m),
        abs(coord.dec.dms.s),
    )


def get_hms_dms_sober(coord: SkyCoord):
    return "{:2.0f} {:02.0f} {:02.2f}  {:2.0f} {:02.0f} {:02.2f}".format(
        coord.ra.hms.h,
        abs(coord.ra.hms.m),
        abs(coord.ra.hms.s),
        coord.dec.dms.d,
        abs(coord.dec.dms.m),
        abs(coord.dec.dms.s),
    )


def get_hms_dms_matplotlib(coord: SkyCoord):
    return r"{:2.0f}$^h$ {:02.0f}$^m$ {:02.2f}$^s$ | {:2.0f}$\degree$ {:02.0f}$'$ {:02.2f}$''$".format(
        coord.ra.hms.h,
        abs(coord.ra.hms.m),
        abs(coord.ra.hms.s),
        coord.dec.dms.d,
        abs(coord.dec.dms.m),
        abs(coord.dec.dms.s),
    )


def get_lesve_coords(coord: SkyCoord):
    return "{:2.0f} {:02.0f} {:02.2f} {:2.0f} {:02.0f} {:02.2f}".format(
        coord.ra.hms.h,
        abs(coord.ra.hms.m),
        abs(coord.ra.hms.s),
        coord.dec.dms.d,
        abs(coord.dec.dms.m),
        abs(coord.dec.dms.s),
    )


def get_pool(processes=cpu_count() - 1, maxtasksperchild=10):
    return mp.Pool(processes, maxtasksperchild=maxtasksperchild)


def add_metadata(stars: List[star_description.StarDescription], metadata: StarMetaData):
    """
    Add a static StarMetaData (or children) object to all stars in the list
    :param stars:
    :param metadata:
    :return:
    """
    for star in stars:
        star.metadata = metadata


def get_stars_with_metadata(
    stars: List[star_description.StarMetaData],
    catalog_name: str,
    exclude: List[str] = [],
) -> List[star_description.StarDescription]:
    # gets all stars which have a catalog of name catalog_name
    assert isinstance(exclude, list) and isinstance(stars, list)
    return list(
        filter(
            partial(metadata_filter, catalog_name=catalog_name, exclude=exclude), stars
        )
    )


def concat_sd_lists(*star_descriptions):
    result = []
    id_set = set()
    for sd_list in star_descriptions:
        assert isinstance(sd_list, list)
        for sd in sd_list:
            if sd.local_id not in id_set:
                result.append(sd)
                id_set.add(sd.local_id)
    return result


# Does this star have a catalog with catalog_name? Used in combination with filter()
def metadata_filter(star: StarDescription, catalog_name, exclude=[]):
    catalogs = star.get_metadata_list()
    return catalog_name in catalogs and len([x for x in exclude if x in catalogs]) == 0


class MetadataSorter:
    pattern = re.compile(r".*?(\d+)$")  # finding the number in our name

    def get_mixed_sort_value(
        self, startuple: Tuple[int, StarDescription], names: List[str]
    ):
        """ gets the value to sort, works with mixed int/str types """
        idx, _ = startuple
        name = names[idx]
        if name is None:
            result = -1, -1
        elif isinstance(name, int) or name.isdigit():
            result = 0, int(name)
        else:
            result = 1, self.get_string_number_part_or_default(name)
        return result

    def get_string_number_part_or_default(
        self, star_name: str, default_value: int = sys.maxsize
    ):
        match = re.match(self.pattern, star_name)
        return int(match.group(1)) if match is not None else default_value

    @staticmethod
    def get_metadata_from_star(
        star: StarDescription, metadata_id: str, warnings: bool = False
    ):
        result = star.get_metadata(metadata_id)
        if result is None and warnings:
            logging.warning(
                f"The metadata {metadata_id} for star {star.local_id} does not exist"
            )
        return result

    @staticmethod
    def get_name_from_metadata(obj, name_var: str, warning: bool = False):
        try:
            return getattr(obj, name_var)
        except AttributeError:
            if warning:
                logging.warning(
                    f"The metadata {obj} does not have a name variable called {name_var}"
                )
            return None

    def __call__(
        self,
        stars: List[StarDescription],
        metadata_id="SITE",
        name_variable="name",
        warnings=True,
    ):
        metadata = [
            MetadataSorter.get_metadata_from_star(x, metadata_id, warnings)
            for x in stars
        ]
        names = [
            MetadataSorter.get_name_from_metadata(x, name_variable, warnings)
            for x in metadata
        ]
        sorted_stars = [
            x[1]
            for x in sorted(
                enumerate(stars), key=partial(self.get_mixed_sort_value, names=names)
            )
        ]
        return sorted_stars


metadata_sorter = MetadataSorter()


def sort_selected(stars: List[StarDescription]) -> List[StarDescription]:
    non_vsx = get_stars_with_metadata(stars, "SITE", exclude=["VSX"])
    # logging.info(f"Non vsx stars: {[x for x in non_vsx]}")
    vsx = get_stars_with_metadata(stars, "VSX")
    # logging.info(f"Vsx stars: {[x for x in non_vsx]}")
    assert len(stars) == len(non_vsx) + len(vsx)
    non_vsx_sorted_stars = metadata_sorter(
        non_vsx, metadata_id="SITE", name_variable="our_name"
    )
    # logging.info(f"Non vsx stars sorted: {[x for x in non_vsx_sorted_stars]}")
    vsx_sorted_stars = metadata_sorter(
        vsx, metadata_id="SITE", name_variable="our_name"
    )
    # logging.info(f"Vsx stars sorted: {[x for x in vsx_sorted_stars]}")
    return non_vsx_sorted_stars + vsx_sorted_stars


def add_star_lists(list1: List[StarDescription], list2: List[StarDescription]):
    ids = [x.local_id for x in list1]
    list2_filtered = [x for x in list2 if x.local_id not in ids]
    return list1 + list2_filtered


def reject_outliers_iqr(df, column, cut=5):
    q1, q3 = np.percentile(df[column], [cut, 100 - cut])
    iqr = q3 - q1
    lower_bound = q1 - (iqr * 1.5)
    upper_bound = q3 + (iqr * 1.5)
    logging.debug(f"q1 {q1} q3 {q3} iqr {iqr} lower {lower_bound} upper {upper_bound}")
    return df[(df[column] < upper_bound) & (df[column] > lower_bound)]


def get_star_or_catalog_name(star: StarDescription, suffix: str = "") -> StarUI:
    extradata = None
    if star.has_metadata("VSX"):
        catalog = star.get_metadata("VSX")
        catalog_name, separation = catalog.name, catalog.separation
        extradata = catalog.extradata
    elif star.has_metadata("SITE"):
        catalog = star.get_metadata("SITE")
        catalog_name, separation = catalog.our_name, catalog.separation
    else:
        catalog_name, separation = star.local_id, None
    filename_no_suff_no_ext = (
        f"{int(catalog_name):05}"
        if isinstance(catalog_name, int) or catalog_name.isdigit()
        else f"{catalog_name}"
    )
    filename_no_ext = f"{filename_no_suff_no_ext}{suffix}"

    filename_orig_no_ext = filename_no_ext
    filename_no_ext = replace_spaces(replace_dots(filename_orig_no_ext))

    filename_orig_no_suff_no_ext = filename_no_suff_no_ext
    filename_no_suff_no_ext = replace_spaces(replace_dots(filename_orig_no_suff_no_ext))
    return StarUI(
        catalog_name,
        separation,
        extradata,
        filename_orig_no_ext,
        filename_no_ext,
        filename_orig_no_suff_no_ext,
        filename_no_suff_no_ext,
    )


def get_star_names(star: StarDescription) -> List[str]:
    def unique_append(alist, new):
        if new not in alist:
            alist.append(new)

    names = []
    if star.has_metadata("VSX"):
        unique_append(names, star.get_metadata("VSX").name)
    if star.has_metadata("SITE"):
        unique_append(names, star.get_metadata("SITE").our_name)
    return names if len(names) > 0 else None


def get_pretty_ucac4_of_sd(star: StarDescription):
    catdata: CatalogData = get_ucac4_of_sd(star)
    return catdata.catalog_id if catdata is not None else "Unknown"


def get_ucac4_of_sd(star: StarDescription) -> CatalogData:
    catdata: CatalogData = star.get_metadata("UCAC4")
    return catdata


def get_full_ucac4_id(ucac4_input: str) -> str:
    """ Takes a partial ucac4 id and makes it complete. E.g.: '233-155284' ==> 'UCAC4 233-155284'"""
    # UCAC4 233-155284
    if len(ucac4_input) >= 10:
        return f"UCAC4 {ucac4_input[-10:]}"
    return None


# replace spaces with dashes
def replace_dots(a_string: str):
    return a_string.replace(".", "-")


# replace spaces with underscores
def replace_spaces(a_string: str):
    return a_string.replace(" ", "_")


# replace spaces with underscores
def replace_underscores(a_string: str):
    return a_string.replace("_", " ")

# if var type is L or period is -1/None then it's not periodic
def is_var_type_aperiodic(var_type, period: float):
    if var_type in ['L', 'None'] or period == -1:
        return True
    return False

# not used


def find_index_of_file(the_dir, the_file, the_filter="*"):
    the_dir = glob.glob(the_dir + "*" + the_filter)
    the_dir.sort()
    indices = [i for i, elem in enumerate(the_dir) if the_file in elem]
    return indices[0]


# not used


def find_file_for_index(the_dir, index, the_filter="*"):
    the_dir = glob.glob(the_dir + the_filter)
    the_dir.sort()
    return the_dir[index]

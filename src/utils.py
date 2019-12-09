import glob
from functools import partial
from os import listdir
from os.path import isfile, join
from subprocess import call

import numpy as np
import star_description
from star_description import StarDescription, StarMetaData
from typing import List, Dict
import multiprocessing as mp
from multiprocessing import cpu_count
import re
import logging


def find_index_of_file(the_dir, the_file, the_filter='*'):
    the_dir = glob.glob(the_dir + "*" + the_filter)
    the_dir.sort()
    indices = [i for i, elem in enumerate(the_dir) if the_file in elem]
    return indices[0]


def find_file_for_index(the_dir, index, the_filter='*'):
    the_dir = glob.glob(the_dir + the_filter)
    the_dir.sort()
    return the_dir[index]


def add_trailing_slash(the_path):
    return join(the_path, '')


def get_files_in_dir(mypath):
    return [f for f in listdir(mypath) if isfile(join(mypath, f))]


# out012345.dat -> 12345
def get_starid_from_outfile(outfile) -> int:
    m = re.search('out(.*).dat', outfile)
    return int(m.group(1).lstrip('0'))


# returns a dict with the local_id as key
def get_star_description_cache(stars: List[StarDescription]) -> Dict[int, StarDescription]:
    cachedict = {}
    for sd in stars:
        cachedict[sd.local_id] = sd
    return cachedict


# filter a list of star descriptions on the presence of a catalog
def catalog_filter(star: StarDescription, catalog_name):
    return star.has_metadata(catalog_name)


def get_hms_dms(coord):
    return "{:2.0f}h {:02.0f}m {:02.2f}s  {:2.0f}d {:02.0f}' {:02.2f}\"" \
        .format(coord.ra.hms.h, abs(coord.ra.hms.m), abs(coord.ra.hms.s),
                coord.dec.dms.d, abs(coord.dec.dms.m), abs(coord.dec.dms.s))


def get_hms_dms_matplotlib(coord):
    return "{:2.0f}$^h$ {:02.0f}$^m$ {:02.2f}$^s$ | {:2.0f}$\degree$ {:02.0f}$'$ {:02.2f}$''$" \
        .format(coord.ra.hms.h, abs(coord.ra.hms.m), abs(coord.ra.hms.s),
                coord.dec.dms.d, abs(coord.dec.dms.m), abs(coord.dec.dms.s))


def get_lesve_coords(coord):
    return "{:2.0f} {:02.0f} {:02.2f} {:2.0f} {:02.0f} {:02.2f}" \
        .format(coord.ra.hms.h, abs(coord.ra.hms.m), abs(coord.ra.hms.s),
                coord.dec.dms.d, abs(coord.dec.dms.m), abs(coord.dec.dms.s))


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


def get_stars_with_metadata(stars: List[star_description.StarMetaData], catalog_name: str,
                            exclude=[]) -> List[star_description.StarDescription]:
    # gets all stars which have a catalog of name catalog_name
    assert isinstance(exclude, list) and isinstance(stars, list)
    return list(filter(partial(metadata_filter, catalog_name=catalog_name, exclude=exclude), stars))


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
    pattern = re.compile(r'.*?(\d+)$')  # finding the number in our name


    def sort_metadata_name(self, stars: List[StarDescription], metadata_id, name_extract, warnings):
        def get_sort_value(star: StarDescription):
            metadata_entry = star.get_metadata(metadata_id)
            number_part = self.get_metadata_name_number_part(
                name_extract(metadata_entry)) if metadata_entry is not None else None
            if warnings and (metadata_entry is None or number_part is None):
                logging.warning(
                    f"Lookup with {metadata_id} id gave name "
                    f"'{name_extract(metadata_entry) if metadata_entry is not None else 'None'}' "
                    f"can't be parsed for sorting, won't be sorted, star: {star}")
                return 0
            return number_part


        return sorted(stars, key=get_sort_value)


    def get_metadata_name_number_part(self, star_name: str):
        match = re.match(self.pattern, star_name)
        return int(match.group(1)) if match is not None else None


    def __call__(self, stars: List[StarDescription], metadata_id='STARFILE', name_extract=lambda x: x.our_name,
                 warnings=True):
        return self.sort_metadata_name(stars, metadata_id, name_extract, warnings)


metadata_sorter = MetadataSorter()


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


def get_star_or_catalog_name(star: StarDescription, suffix: str):
    extradata = None
    if star.has_metadata("VSX"):
        catalog = star.get_metadata("VSX")
        catalog_name, separation = catalog.name, catalog.separation
        extradata = catalog.extradata
    elif star.has_metadata("OWNCATALOG"):
        catalog = star.get_metadata("OWNCATALOG")
        catalog_name, separation = catalog.name, catalog.separation
    else:
        catalog_name, separation = None, None
    filename_no_ext = f"{catalog_name}{suffix}" if catalog_name is not None else f"{star.local_id:05}{suffix}"
    return catalog_name, separation, extradata, replace_spaces(filename_no_ext)


# replace spaces with underscores
def replace_spaces(a_string: str):
    return a_string.replace(' ', '_')

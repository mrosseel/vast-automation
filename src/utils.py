import glob
from os import listdir
from os.path import isfile, join
from star_description import StarDescription
from typing import List, Dict

def find_index_of_file(the_dir, the_file, the_filter='*'):
    the_dir = glob.glob(the_dir + "*"+the_filter)
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


# returns a dict with the local_id as key
def get_star_description_cache(stars: List[StarDescription]) -> Dict[int, StarDescription]:
    cachedict = {}
    for sd in stars:
        cachedict[sd.local_id] = sd
    return cachedict

# filter a list of star descriptions on the presence of a catalog
def catalog_filter(star: StarDescription, catalog_name):
    return star.has_catalog(catalog_name)

def get_hms_dms(coord):
    return "{:2.0f}h {:02.0f}m {:02.2f}s | {:2.0f}d {:02.0f}' {:02.2f}\"" \
        .format(coord.ra.hms.h, abs(coord.ra.hms.m), abs(coord.ra.hms.s),
                coord.dec.dms.d, abs(coord.dec.dms.m), abs(coord.dec.dms.s))

def get_lesve_coords(coord):
    return "{:2.0f} {:02.0f} {:02.2f} {:2.0f} {:02.0f} {:02.2f}" \
        .format(coord.ra.hms.h, abs(coord.ra.hms.m), abs(coord.ra.hms.s),
                coord.dec.dms.d, abs(coord.dec.dms.m), abs(coord.dec.dms.s))

# Reads all photometry files (pht) and returns a big numpy array with all relevant data in it

# JD V-C s1 V s2 C s3
# Aperture: 2, Filter: V
# 2457657.5088310 -0.50728 0.10291 16.65794 0.05604 17.16522 0.08631

from init_loader import init
from reading import trash_and_recreate_dir
from tqdm import tqdm
from functools import partial
import multiprocessing as mp
from read_pht import read_pht_file
import glob
import numpy as np
from typing import List
import logging

Vector = List[float]
STAR_DATA_MB = 1.5274047851562502e-05 # the size of the data of one star

# Returns: jd[nrfiles], fwhm[nrfiles], star_result[nrfiles, nrstars, 2]
#   the func is now passed a range of stars to process. We now need to iterate through all the pht files
#   and extract info about ONLY those stars, then write them to file. Maybe dense numpy lists with the range as index?
def read_photometry(star_list_1, apertureidx: int, matched_files=None, fake_references=False):
    if matched_files is None:
        matched_files = glob.glob(init.matchedphotometrydir+"*.pht") # todo extract dir, pattern
    nrfiles = len(matched_files)
    nrstars = len(star_list_1)
    star_range_0 = np.array(star_list_1) - 1
    star_result_ = np.full([nrfiles, nrstars, 2],np.inf, dtype=float)
    fwhm = np.full([nrfiles, 3], np.inf, dtype=float)
    jd = np.full([nrfiles], np.inf, dtype=float)
    nrstars = np.full([nrfiles], np.inf, dtype=int)
    nr_threads=init.nr_threads*2
    pbar = tqdm(total=len(matched_files))
    pool = mp.Pool(nr_threads, maxtasksperchild=None)
    func = partial(read_pht, star_range_0=star_range_0, apertureidx=apertureidx, fake_references=fake_references)

    for fileidx, jd_, fwhm_, nrstars_, collected in pool.imap_unordered(func, enumerate(matched_files)):
        jd[fileidx] = jd_
        fwhm[fileidx] = fwhm_
        nrstars[fileidx] = nrstars_
        star_result_[fileidx] = collected
        pbar.update(1)
    pbar.close()
    pool.close()
    pool.join()
    return jd, fwhm, nrstars, star_result_

def read_pht(matched_files_tuple, star_range_0, apertureidx: int, fake_references=False):
    fileidx = matched_files_tuple[0]
    file_entry = matched_files_tuple[1]
    # open the file for reading
    with open(file_entry, mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        photheader, _, nrstars, stars, stardata = read_pht_file(file_entry, fileContent, only_apertureidx=int(apertureidx))
        collect = np.full([len(init.star_list), 2],np.inf, dtype=float) # we don't use nrstars because this could be e.g. 1000, but with stars which have id's > 1000

        fwhm = [photheader.fwhm_exp, photheader.fwhm_mean, photheader.fwhm_err]
        jd = photheader.jd

        # for every star
        for staridx, starentry in enumerate(stardata):
            # if stars[staridx].ref_id-1 == 13:
            #   print("readandcollectpht:", [starentry.mag, starentry.err], stars[staridx].ref_id-1)
            ref_id_0 = stars[staridx].ref_id - 1 # if photometry then all is -1

            if ref_id_0 is -2: # -1 - 1
                if fake_references:
                    ref_id_0 = staridx
                else:
                    continue
            try:
                collect[ref_id_0] = [starentry.mag, starentry.err]
            except:
                logging.error(f"staridx: {staridx}, ref_id_0: {ref_id_0}, shape: {collect.shape}")
        collected = collect[star_range_0]
        return fileidx, jd, fwhm, nrstars, collected

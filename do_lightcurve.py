# JD V-C s1 V s2 C s3
# Aperture: 2, Filter: V
# 2457657.5088310 -0.50728 0.10291 16.65794 0.05604 17.16522 0.08631

import init
from reading import trash_and_recreate_dir
from reading import reduce_star_list
import tqdm
from functools import partial
from subprocess import call
import multiprocessing as mp
from multiprocessing.pool import ThreadPool
from do_photometry import read_pht_file
import glob
import numpy as np
import math
import sys
from typing import List

photometry=None
preamble=None
Vector = List[float]
MAX_MAG = 99.99999
MAX_ERR = 9.99999
STAR_DATA_MB=1.5274047851562502e-05 # the size of the data of one star
star_result = None

# star_list is ignored for the moment or always if it's fast enough
def main(star_list_1, check_stars_1, aperture, apertureidx, is_resume):
    if not is_resume:
        trash_and_recreate_dir(init.lightcurvedir)
    else:
        star_list_1 = reduce_star_list(star_list_1, init.lightcurvedir)

    global preamble
    preamble = init_preamble(aperture, check_stars_1)
    #print("Preamble is", preamble)
    # chop up the stars into parts of CHUNK_SIZE each
    matched_files = glob.glob(init.matchedphotometrydir+"*.pht") # todo extract dir, pattern

    # calculate possible batch size
    chunk_size = min((init.free_memory_GB / (len(star_list_1)*len(matched_files)*STAR_DATA_MB/1024))*len(star_list_1), len(star_list_1))
    print("Chunk size is", chunk_size)
    star_ranges = chunks(star_list_1, chunk_size)

    # pool = mp.Pool(init.nr_threads, maxtasksperchild=None)
    pool = ThreadPool(1)
    func = partial(read_and_write_star_data, check_stars_1=check_stars_1, aperture=aperture, apertureidx=apertureidx, matched_files=matched_files)
    #print("Writing star lightcurves for", len(star_list_1), "stars into", init.lightcurvedir)
    pool.imap_unordered(func, star_ranges)

#   the func is now passed a range of stars to process. We now need to iterate through all the pht files
#   and extract info about ONLY those stars, then write them to file. Maybe dense numpy lists with the range as index?
def read_and_write_star_data(star_range_1, check_stars_1, aperture, apertureidx, matched_files):
    #print(f"Read and write star data with range {star_range_1}")
    jd, fwhm, star_result_ = read_star_data(star_range_1, matched_files, apertureidx)
    global star_result
    star_result = star_result_
    pool = mp.Pool(init.nr_threads*2, maxtasksperchild=None)
    func = partial(write_lightcurve, check_stars_1=check_stars_1, aperture=aperture, apertureidx=apertureidx, jd=jd, fwhm=fwhm)

    #print("Writing star lightcurves for", len(star_list_1), "stars into", init.lightcurvedir)
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_range_1), total=len(star_range_1), desc='Writing lightcurve'):
        pass

def read_star_data(star_range_1, matched_files, apertureidx):
    #print(f"Read star data with range {star_range_1}")
    nrfiles = len(matched_files)
    nrstars = len(star_range_1)
    star_range_0 = np.array(star_range_1) - 1
    star_result = np.empty([nrfiles, nrstars, 2],dtype=float)
    fwhm = np.empty([nrfiles, 3], dtype=float)
    jd = np.empty([nrfiles], dtype=float)

    # gather all the data
    for fileidx, file_entry in tqdm.tqdm(enumerate(matched_files), total=len(matched_files), desc='Read pht'):
        fileContent = None
        # open the file for reading
        with open(file_entry, mode='rb') as file: # b is important -> binary
            fileContent = file.read()
            jd_, fwhm_, collect = read_and_collect_pht(file_entry, fileContent, apertureidx)
            jd[fileidx] = jd_
            fwhm[fileidx] = fwhm_
            # print("Collect shape", collect.shape)
            print(f"starrangezero len: {len(star_range_0)}, collect[] len: {len(collect[apertureidx])}")

            collected = collect[apertureidx][star_range_0]
            # print("Collected is:", collected, collected.shape)
            star_result[fileidx] = collected
            # print("star result is:", star_result.shape)
            # print("star[fil]:", star_result[fileidx])
            # print("star[fil][0]:", star_result[fileidx][0])
            # print("star[fil][0][0]:", star_result[fileidx][0][0])
    #print("star[0][0][0]:", star_result[0][0])
    #print("read star data nbytes:", star_result.nbytes, jd.nbytes, fwhm.nbytes)
    return jd, fwhm, star_result

def read_and_collect_pht(file_entry, fileContent, apertureidx):
    # print(f"Read and collect pht: {file_entry}")
    # Star = namedtuple('Star', 'id, ref_id, x, y, skymed, skysig, fwhm')
    # StarData = namedtuple('StarData', 'mag, err, code')
    photheader, _, nrstars, stars, stardata = read_pht_file(file_entry, fileContent, only_apertureidx=apertureidx)
    # nr of apertures = len(stardata[0]
    nrapertures = len(stardata[0])
    collect = np.empty([nrapertures, nrstars, 2],dtype=float)

    # print("\tDate from header:",photheader.jd, "fwhm:", photheader.fwhm_mean)
    fwhm = [photheader.fwhm_exp, photheader.fwhm_mean, photheader.fwhm_err]
    jd = photheader.jd
    # logging.debug(f"the result is {result[0]}")

    # for every star
    for staridx, starentry in enumerate(stars):
        # for every aperture
        for apidx, aperturedata in enumerate(stardata[staridx]):
            if apidx != apertureidx:
                continue
            if aperturedata.mag == 0:
                print(aperturedata)
            collect[apidx][starentry.ref_id-1] = [aperturedata.mag, aperturedata.err]
    return jd, fwhm, collect

def write_lightcurve(star_1: int, check_stars_1: Vector, aperture: float, apertureidx: int, jd: float, fwhm: float):
    #print(f"Write lightcurve for star:", star_1)
    check_stars = join_check_stars(check_stars_1, star_1)
    all_stars_1 = [star_1]
    all_stars_1 = all_stars_1 + check_stars
    all_stars_0 = np.array(all_stars_1) - 1
    #print("write lightcurve all_stars_0", all_stars_0)
    # print("check stars:", check_stars, "allstars:", all_stars_1)
    # print("photometry:", photometry)
    lines = []
    lines.append(preamble)
    #print("nr of files:", nrfiles)
    # for every file
    sorted_jd = np.argsort(jd) # argsort the julian date so lines are inserted in the correct order
    for fileidx in sorted_jd:
        line = ""
        date = f"{jd[fileidx]:.7f}"
        line = line + date
        V = min(MAX_MAG, star_result[fileidx][all_stars_0[0]][0])
        Verr = min(MAX_ERR, star_result[fileidx][all_stars_0[0]][1])
        C = min(MAX_MAG, star_result[fileidx][all_stars_0[1]][0])
        Cerr = min(MAX_ERR, star_result[fileidx][all_stars_0[1]][1])
        line = line + f" {(V-C):.5f} {math.sqrt(Verr**2 + Cerr**2):.5f}"
        for allstaridx, _ in enumerate(all_stars_1):
            line = line + f" {min(MAX_MAG, star_result[fileidx][allstaridx][0]):.5f} {min(MAX_ERR, star_result[fileidx][allstaridx][1]):.5f}"
        lines.append(line)

    with open(init.lightcurvedir + 'curve_' + str(star_1).zfill(5) + ".txt", 'wt') as f:
        f.write('\n'.join(lines)+'\n')

# the static first part of the file
def init_preamble(aperture, check_stars_list):
    preamble = "JD V-C s1 V s2"
    checkstarcount = 1
    count = 3
    for _ in check_stars_list:
        preamble = preamble + f" C{checkstarcount} s{count}"
        checkstarcount = checkstarcount + 1
        count = count + 1
    preamble = preamble + f"\nAperture: {aperture}, Filter: V, Check stars: {check_stars_list}"
    return preamble

def join_check_stars_string(check_stars, exclude_star):
    check_stars = filter(lambda star: star != exclude_star, check_stars)
    check_stars_string = ','.join(map(str, check_stars))
    return check_stars_string

def join_check_stars(check_stars, exclude_star):
    check_stars = list(filter(lambda star: star != exclude_star, check_stars))
    return check_stars

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

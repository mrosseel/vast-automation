# JD V-C s1 V s2 C s3
# Aperture: 2, Filter: V
# 2457657.5088310 -0.50728 0.10291 16.65794 0.05604 17.16522 0.08631

import init
from reading import trash_and_recreate_dir
from reading import reduce_star_list
from tqdm import tqdm
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
def main(star_list_1_orig, check_stars_1, aperture, apertureidx, is_resume):
    if not is_resume:
        trash_and_recreate_dir(init.lightcurvedir)
        star_list_1 = star_list_1_orig
    else:
        star_list_1 = reduce_star_list(star_list_1_orig, init.lightcurvedir)

    global preamble
    preamble = init_preamble(aperture, check_stars_1)
    #print("Preamble is", preamble)
    # chop up the stars into parts of CHUNK_SIZE each
    matched_files = glob.glob(init.matchedphotometrydir+"*.pht") # todo extract dir, pattern

    # calculate possible batch size
    chunk_one = (init.free_memory_GB / (len(star_list_1)*len(matched_files)*STAR_DATA_MB/1024))*len(star_list_1)
    chunk_two = len(star_list_1)
    chunk_size = min(chunk_one, chunk_two)
    star_ranges = chunks(star_list_1, chunk_size)
    print(f"Reading {chunk_size}/{len(star_list_1)} stars per batch.")

    # pool = mp.Pool(init.nr_threads, maxtasksperchild=None)
    pool = ThreadPool(1)
    func = partial(read_and_write_star_data, check_stars_1=check_stars_1, aperture=aperture, apertureidx=apertureidx, matched_files=matched_files)
    #print("Writing star lightcurves for", len(star_list_1), "stars into", init.lightcurvedir)
    for _ in pool.imap_unordered(func, star_ranges):
        pass

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
    for _ in tqdm(pool.imap_unordered(func, star_range_1), total=len(star_range_1), desc='Writing lightcurve'):
        pass

def read_star_data(star_range_1, matched_files, apertureidx):
    #print(f"Read star data with range {star_range_1}")
    nrfiles = len(matched_files)
    nrstars = len(star_range_1)
    star_range_0 = np.array(star_range_1) - 1
    star_result_ = np.empty([nrfiles, nrstars, 2],dtype=float)
    fwhm = np.empty([nrfiles, 3], dtype=float)
    jd = np.empty([nrfiles], dtype=float)

    pbar = tqdm(total=len(matched_files))
    pool = mp.Pool(init.nr_threads*2, maxtasksperchild=None)
    func = partial(read_pht, star_range_0=star_range_0, apertureidx=apertureidx)

    for fileidx, jd_, fwhm_, collected in pool.imap_unordered(func, enumerate(matched_files)):
        jd[fileidx] = jd_
        fwhm[fileidx] = fwhm_
        star_result_[fileidx] = collected
        pbar.update(1)
    pbar.close()
    return jd, fwhm, star_result_

def read_pht(matched_files_tuple, star_range_0, apertureidx):
    fileContent = None
    fileidx = matched_files_tuple[0]
    file_entry = matched_files_tuple[1]
    # open the file for reading
    with open(file_entry, mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        jd_, fwhm_, collect = read_and_collect_pht(file_entry, fileContent, apertureidx, len(star_range_0))
        # print("Collect shape", collect.shape)
        # print(f"fileidx: {fileidx}, star_range_0 len: {len(star_range_0)}, collect[apixd] len: {len(collect[apertureidx])}, {star_range_0}")
        # print("collect apertureidx", apertureidx, collect[apertureidx][0:14])
        collected = collect[star_range_0]
        # print("collected:", collected[13])
        # print("Collected is:", collected, collected.shape)
        # star_result[fileidx] = collected
        return fileidx, jd_, fwhm_, collected

def read_and_collect_pht(file_entry, fileContent, apertureidx:int, star_range_length):
    # print(f"Read and collect pht: {file_entry}")
    # Star = namedtuple('Star', 'id, ref_id, x, y, skymed, skysig, fwhm')
    # StarData = namedtuple('StarData', 'mag, err, code')
    photheader, _, _, stars, stardata = read_pht_file(file_entry, fileContent, only_apertureidx=apertureidx)
    collect = np.empty([star_range_length, 2],dtype=float) # we don't use nrstars because this could be e.g. 1000, but with stars which have id's > 1000

    # print("\tDate from header:",photheader.jd, "fwhm:", photheader.fwhm_mean)
    fwhm = [photheader.fwhm_exp, photheader.fwhm_mean, photheader.fwhm_err]
    jd = photheader.jd
    # logging.debug(f"the result is {result[0]}")

    # for every star
    for staridx, starentry in enumerate(stardata):
        # if stars[staridx].ref_id-1 == 13:
        #   print("readandcollectpht:", [starentry.mag, starentry.err], stars[staridx].ref_id-1)
        collect[stars[staridx].ref_id-1] = [starentry.mag, starentry.err]
    return jd, fwhm, collect

def write_lightcurve(star_1: int, check_stars_1: Vector, aperture: float, apertureidx: int, jd: float, fwhm: float):
    #print(f"Write lightcurve for star:", star_1)
    check_stars_1 = join_check_stars(check_stars_1, star_1)
    check_stars_0 = np.array(check_stars_1) - 1
    star_0 = star_1 - 1

    lines = [preamble]
    #print("nr of files:", nrfiles)
    # for every file
    sorted_jd = np.argsort(jd) # argsort the julian date so lines are inserted in the correct order
    for fileidx in sorted_jd:
        line = f"{jd[fileidx]:.7f}" # start the line with the julian date
        V = star_result[fileidx][star_0][0]
        Verr = min(MAX_ERR, star_result[fileidx][star_0][1])
        if not is_valid(V, Verr): continue
        # if V < 1 or V > 99:
        #     print(f"Strange V for fileidx {fileidx}, V:{V}, Verr: {Verr}, Star_1: {star_1}")
        C, Cerr = calculate_synthetic_c(star_result[fileidx], check_stars_0)
        if C == 0 and Cerr == 0: continue # abort if one of the comparison stars is not available
        Cerr = min(MAX_ERR, Cerr)
        # if star_1 == 14:
        #     print(f"fileidx {fileidx}, V:{V}, V2:{star_result[fileidx][star_0][0]} Verr: {Verr}, Star_1: {star_1}, C: {C}, Cerr: {Cerr}")

        linedata = [(V-C, math.sqrt(Verr**2 + Cerr**2)), (V, Verr), (C, Cerr)] \
            + [(star_result[fileidx][checkstar_0][0], star_result[fileidx][checkstar_0][1]) for checkstar_0 in check_stars_0]
        # print(linedata)
        for tuple in linedata:
            line += f" {min(MAX_MAG, tuple[0]):.5f} {min(MAX_ERR, tuple[1]):.5f}"
        lines.append(line)

    with open(init.lightcurvedir + 'curve_' + str(star_1).zfill(5) + ".txt", 'wt') as f:
        # for l in lines: f.write('%s\n' % l)
        f.write('\n'.join(lines)+'\n')

# the static first part of the file
def init_preamble(aperture, check_stars_list):
    preamble = "JD V-C s1 V s2 C s3"
    checkstarcount = 1
    count = 4
    for _ in check_stars_list:
        preamble = preamble + f" C{checkstarcount} s{count}"
        checkstarcount = checkstarcount + 1
        count = count + 1
    preamble = preamble + f"\nAperture: {aperture}, Filter: V, Check stars: {check_stars_list}"
    return preamble

# /* Comparison star */
# if (lc->comp.count==1) {
# 	if (comp[0].valid) {
# 		cmag = comp[0].mag;
# 		cerr = comp[0].err;
# 		comp_ok = 1;
# 	} else {
# 		cerr = cmag = 0.0;
# 		comp_ok = 0;
# 	}
# } else {
# 	cmag = cerr = 0.0;
# 	n = 0;
# 	for (i=0; i<lc->comp.count; i++) {
# 		if (comp[i].valid) {
# 			cmag += pow(10.0, -0.4*comp[i].mag);
# 			cerr += comp[i].err;
# 			n++;
# 		}
# 	}
# 	if (n==lc->comp.count) {
# 		cmag = -2.5*log10(cmag/n);
# 		cerr = (cerr/n)/sqrt((double)n);
# 		comp_ok = 1;
# 	} else {
# 		cerr = cmag = 0.0;
# 		comp_ok = 0;
# 	}
# }
def calculate_synthetic_c(star_result_file, check_stars_0):
    if len(check_stars_0) == 1:
        cmag = star_result_file[check_stars_0[0]][0]
        cerr = star_result_file[check_stars_0[0]][1]
        if is_valid(cmag, cerr):
            return cmag, cerr
        else:
            return 0,0

    cmag, cerr, valid = 0, 0, 0
    nrstars = len(check_stars_0)
    for entry in check_stars_0:
        mag = star_result_file[entry][0]
        err = star_result_file[entry][1]
        if is_valid(mag, err):
            cmag += math.pow(10, -0.4*mag)
            cerr += err
            valid += 1
    if valid == nrstars:
        cmag = -2.5*math.log10(cmag/nrstars)
        cerr = (cerr/nrstars)/math.sqrt(nrstars)
    else:
        cmag, cerr = 0, 0
    return cmag, cerr

def is_valid(mag, err):
    return not np.isnan(mag) and not np.isnan(err) and mag < MAX_MAG and err < MAX_ERR

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

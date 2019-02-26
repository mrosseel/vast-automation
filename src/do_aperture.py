import glob
import sys
from tqdm import tqdm
import numpy as np
import logging
import scipy
import init
from read_pht import read_pht_file

# Select files conforming to the match_pattern using percentage which is between 0 and 1
def select_files(the_dir, match_pattern, percentage=1):
    matched_files = glob.glob(the_dir+match_pattern)
    desired_length = max(1, int(len(matched_files) * percentage))
    np.random.seed(42) # for the same percentage, we always get the same selection
    selected_files = np.random.choice(matched_files, size=desired_length, replace=False).tolist()
    return selected_files

def convert_to_aperture_only(apertures):
    return apertures[1::2]

def main(the_dir, match_files='match*.pht', percentage=0.1):
    files = select_files(the_dir, match_files, percentage)
    nrfiles = len(files)
    detectablestars = len(init.star_list) # this is the maximum of star id's, regardless of how many are detected on an image
    logging.debug(files)
    # pre process
    apertures = []
    with open(files[0], mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        photheader, apertures, nrstars , _, _ = read_pht_file(files[0], fileContent)
        apertures = convert_to_aperture_only(apertures)
    file.close()
    logging.info(f"Apertures: {apertures}, nr of files: {nrfiles}")
    collect = np.empty([len(apertures), detectablestars, nrfiles, 2],dtype=float)
    fwhm = np.empty([nrfiles, 3], dtype=float)
    jd = np.empty([nrfiles], dtype=float)
    # for all files
    for fileidx, entry in enumerate(files):
        fileContent = None
        # open the file for reading
        with open(entry, mode='rb') as file: # b is important -> binary
            fileContent = file.read()
            print(f"{fileidx}/{nrfiles}: {entry}")
            photheader, apertures, nrstars, stars, stardata = read_pht_file(entry, fileContent)
            apertures = convert_to_aperture_only(apertures)
            print("\tDate from header:",photheader.jd, "fwhm:", photheader.fwhm_mean)
            fwhm[fileidx] = [photheader.fwhm_exp, photheader.fwhm_mean, photheader.fwhm_err]
            jd[fileidx] = photheader.jd
            # logging.debug(f"the result is {result[0]}")

            # for every star
            for staridx, starentry in enumerate(stars):
                # for every aperture
                for apidx, aperturedata in enumerate(stardata[staridx]):
                    if aperturedata.mag == 0:
                        print(aperturedata)
                    collect[apidx][starentry.ref_id-1][fileidx] = [aperturedata.mag, aperturedata.err]
                    # collect[apidx][starentry.ref_id-1][fileidx][1] = aperturedata.err
                    if collect[apidx][starentry.ref_id-1][fileidx][0] == 0:
                        print("nul")

    print("Nr of stars for apertureidx 2:", len(collect[2]))
    print("Mb used:", collect.nbytes/1024/1024)
    print("Entry for First aperture, first star:", collect[0][0])
    print("Shape for collect:", collect.shape)
    print(scipy.stats.mstats.describe(np.ma.masked_invalid(collect[2]), axis=0))
    stddevs = np.empty([len(apertures), detectablestars])
    counts = np.empty([len(apertures), detectablestars])
    import warnings
    warnings.simplefilter('error', UserWarning)
    logging.info("Calculating stddevs...")
    # pool = mp.Pool(init.nr_threads*2, maxtasksperchild=None)
    # func = partial(write_lightcurve, check_stars_1=check_stars_1, aperture=aperture, apertureidx=apertureidx, jd=jd, fwhm=fwhm)
    # logging.debug("Writing star lightcurves for", len(star_list_1), "stars into", init.lightcurvedir)
    # for _ in tqdm(pool.imap_unordered(func, star_list_1), total=len(star_list_1), desc='Writing lightcurve'):
    #     pass

    for apidx in tqdm(range(len(collect)), desc="Calculating stddevs"):
        for staridx in range(len(collect[apidx])):
            # print(apidx, staridx)
            masked_collect = np.ma.masked_invalid(collect[apidx][staridx])
            # print("Shape of the thing we calculate stddev on:", collect[apidx][staridx].shape, collect[apidx][staridx])
            std = masked_collect.std(axis=0)[0] # take the stddev of the magnitude
            count= masked_collect.count()
            # print("count", count)
            # print("std:", std, type(std), std is np.ma.masked)
            stddevs[apidx, staridx] = std if not std is np.ma.masked else sys.float_info.max
            counts[apidx, staridx] = count

    logging.info(stddevs[0,0:20])
    logging.info(f"Shape for stddevs: {stddevs.shape}")
    logging.info(np.argmin(stddevs[0]))
    # for idx in range(len(stddevs)):
    #     print(stddevs[idx].min(), stddevs[idx].max())
    #     print(idx, stddevs[idx].sum())
    median_erroradded = np.median(np.add(np.take(fwhm, 1, axis=1), np.take(fwhm, 2, axis=1)))
    median_multiply = np.median(np.take(fwhm, 1, axis=1))*1.75

    apertureidx = np.abs(apertures - median_multiply).argmin()
    logging.info(f"FWHM median: {median_multiply} aperture chosen is: {apertures[apertureidx]}")
    # winners = np.argwhere(np.amax(counts[apertureidx])).flatten()
    winners = np.argwhere(counts[apertureidx][counts[apertureidx] > np.max(counts[apertureidx])*0.9]).flatten()
    # winners = counts[apertureidx][counts[apertureidx] > np.max(counts[apertureidx])*0.9]
    logging.info(f"Winners: {winners}")
    winnerslice = stddevs[apertureidx][winners]
    nrtopstars = min(3, len(winnerslice))
    compstars_0 = np.argpartition(winnerslice, range(nrtopstars))[:nrtopstars]

    #compstar_0 = np.argmin(stddevs[apertureidx], axis=0)
    logging.info(f"Compstars_0 with minimum stdev in the chosen aperture: {compstars_0}")
    logging.info(f"Compstars stddev: {stddevs[apertureidx][compstars_0]:.8f}")
    logging.info(f"Compstars counts: {counts[apertureidx][compstars_0]}")
    return stddevs, collect, apertures, apertureidx, fwhm, jd, (compstars_0+1).tolist()

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    read_pht_file(sys.argv[1], open(sys.argv[1],'rb').read())

import glob
import sys
from tqdm import tqdm
import numpy as np
import logging
import scipy.stats
import init
from reading import file_selector
from read_pht import read_pht_file


def convert_to_aperture_only(apertures):
    return apertures[1::2]

def gather_data(match_file_list):
    nrfiles = len(match_file_list)
    detectablestars = len(init.star_list) # this is the maximum of star id's, regardless of how many are detected on an image
    # pre process
    apertures = []
    with open(match_file_list[0], mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        photheader, apertures, nrstars , _, _ = read_pht_file(match_file_list[0], fileContent)
        apertures = convert_to_aperture_only(apertures)
    file.close()
    logging.info(f"Apertures: {apertures}, nr of files: {nrfiles}")
    collect = np.full([len(apertures), detectablestars, nrfiles, 2],np.inf, dtype=float)
    fwhm = np.full([nrfiles, 3], np.inf, dtype=float)
    jd = np.full([nrfiles], np.inf, dtype=float)
    # for all files
    for fileidx, entry in enumerate(match_file_list):
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
    return jd, fwhm, collect, apertures

def process(jd, fwhm, collect, apertures):
    detectablestars = len(init.star_list) # this is the maximum of star id's, regardless of how many are detected on an image
    nrapertures = len(apertures)
    # print(scipy.stats.mstats.describe(np.ma.masked_invalid(collect[2]), axis=0))
    stddevs = np.full([nrapertures, detectablestars], np.inf)
    counts = np.full([nrapertures, detectablestars], np.inf)
    logging.info("Calculating stddevs...")
    # pool = mp.Pool(init.nr_threads*2, maxtasksperchild=None)
    # func = partial(write_lightcurve, check_stars_1=check_stars_1, aperture=aperture, apertureidx=apertureidx, jd=jd, fwhm=fwhm)
    # logging.debug("Writing star lightcurves for", len(star_list_1), "stars into", init.lightcurvedir)
    # for _ in tqdm(pool.imap_unordered(func, star_list_1), total=len(star_list_1), desc='Writing lightcurve'):
    #     pass

    for apidx in tqdm(range(nrapertures), desc="Calculating stddevs"):
        for staridx in range(len(collect[apidx])):
            # print(apidx, staridx)
            masked_collect = np.ma.masked_invalid(collect[apidx][staridx])
            # print("Shape of the thing we calculate stddev on:", collect[apidx][staridx].shape, collect[apidx][staridx])
            std = masked_collect.std(axis=0)[0] # take the stddev of the magnitude
            # TODO this is the reason that counts are off, this counts mag and err together
            count= masked_collect.count()
            # print("count", count)
            # print("std:", std, type(std), std is np.ma.masked)
            stddevs[apidx, staridx] = std if not std is np.ma.masked else np.nan
            counts[apidx, staridx] = count

    logging.info(stddevs[0,0:20])
    logging.info(f"Shape for stddevs: {stddevs.shape}")
    logging.info(np.argmin(stddevs[0]))
    # for idx in range(len(stddevs)):
    #     print(stddevs[idx].min(), stddevs[idx].max())
    #     print(idx, stddevs[idx].sum())

    return jd, fwhm, collect, stddevs, counts

def select_aperture_and_compstars(jd, fwhm, collect, apertures, stddevs, counts):
    apertureidx = select_aperture(fwhm, apertures)
    compstars_1 = select_compstars(apertureidx, stddevs, counts)
    return stddevs, collect, apertures, apertureidx, fwhm, jd, compstars_1

def select_aperture(fwhm, apertures):
    # taking the median fwhm
    median_fwhm = np.median(np.take(fwhm, 1, axis=1))
    # using the lesve heuristic of multiplying the fwhm with 1,75 or 2
    median_multiply = median_fwhm*1.75
    # taking the index of the aperture closest to fwhm * 1.75
    apertureidx = np.abs(apertures - median_multiply).argmin()
    logging.info(f"FWHM median: {median_fwhm}, multiplied: {median_multiply}, aperture chosen is: {apertures[apertureidx]}")
    return apertureidx

def select_compstars(apertureidx, stddevs, counts):
    selectedcounts = counts[apertureidx]
    # winners = np.argwhere(np.amax(selectedcounts)).flatten()
    logging.info(f"Before selection on counts: {len(selectedcounts)}")
    logging.info(f"counts threshold {np.max(selectedcounts)*0.9}")
    count_mask = selectedcounts < np.max(selectedcounts)*0.9
    logging.info(f"Count mask: {count_mask}")
    logging.info(f"After selection on counts: {np.sum(count_mask)}")
    #    logging.info(f"Counts: {sorted(selectedcounts)}")

    masked_std = np.ma.masked_array(stddevs[apertureidx], count_mask)
    logging.info(f"Masked std: {masked_std}")
    logging.info(f"masked std len: {masked_std.count()}")

    nrtopstars = min(10, len(masked_std))
    compstars_0 = masked_std.argsort(fill_value=99999)[:nrtopstars]
    logging.info(f"compstars: {compstars_0}, len: {len(compstars_0)}")

    #compstar_0 = np.argmin(stddevs[apertureidx], axis=0)
    logging.info(f"Compstars_0 with minimum stdev in the chosen aperture: {compstars_0}")
    logging.info(f"Compstars stddev: {stddevs[apertureidx][compstars_0]}")
    logging.info(f"Compstars counts: {selectedcounts[compstars_0]}")
    # for star in compstars_0:
    #     logging.info(f"Error for compstar_0:{star} is \t {errors[star].median()}")
    return (compstars_0+1).tolist()

def main(the_dir, match_files='match*.pht', percentage=0.1):
    jd, fwhm, collect, apertures = gather_data(file_selector(init.matchedphotometrydir, 'match*.pht', init.aperture_find_percentage))
    jd, fwhm, collect, stddevs, counts = process(jd, fwhm, collect, apertures)
    return select_aperture_and_compstars(jd, fwhm, collect, apertures, stddevs, counts)

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    main()


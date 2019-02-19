import os
import struct
from collections import namedtuple
import glob
import init
import sys
import numpy as np
import logging


def read_pht_file(the_file, fileContent, only_apertureidx=-1):
    logging.debug(f"The file:{the_file}")
    MAGIC_STRING = "C-Munipack photometry file\r\n" # length is 28
    logging.debug(f'size of file {len(fileContent)}')
    ############################### HEADER ###############################################
    PhotHeader = namedtuple('PhotHeader', 'magic, revision, headerlen, width, height, jd, filter, exptime, ccdtemp, origin, year, mon, mday, hour, min, sec, \
                    range0, range1, gain, rnoise, fwhm_exp, fwhm_mean, fwhm_err, threshold, sharpness0, sharpness1, roundness0, roundness1, matched, match_rstars, \
                    match_istars, match_mstars, match_clip, offset0, offset1, objectdesig, objectra, objectdec, locationdesignation, \
                    locationlon, locationlat, trafoxx, trafoxy, trafox0, trafoyx, trafoyy, trafoy0 \
                    ')
    headerformat = "<"+str(len(MAGIC_STRING))+"sxxxx4id70sdd70sH5Bxddddddddddddiiii3d70sdd70sdddddddd"
    headerlength = struct.calcsize(headerformat)
    logging.debug(f"header length is {headerlength}")
    unp = struct.unpack(headerformat, fileContent[:headerlength])
    photheader = PhotHeader._make(unp)
    for name, entry in photheader._asdict().items():
        logging.debug(f"{name}: {entry}")
    ############################### WCS ###############################################
    wcsformat= "<i"
    wcssize = struct.calcsize(wcsformat)
    wcsdatalength = struct.unpack("<i", fileContent[headerlength:headerlength+wcssize])[0]
    wcslength = headerlength + wcssize + wcsdatalength
    logging.debug(f"wcs data length: {wcsdatalength}")
    ############################### APERTURES ###############################################
    # typedef struct _CmpackPhtAperture
    # {
    # 	int id;								/**< Aperture identifier */
    # 	double radius;						/**< Radius in pixels */
    # } CmpackPhtAperture;

    apertureformat1 = "<i"
    ap1length = struct.calcsize(apertureformat1)
    start = wcslength
    try:
        napertures = struct.unpack(apertureformat1, fileContent[start:start+ap1length])[0]
    except:
        print("error reading", the_file)
        print("'",fileContent[start:start+ap1length],"'")
    logging.debug(f"napertures: {napertures}")
    apertureformat2 = "<" + napertures * "id"
    ap2length = struct.calcsize(apertureformat2)
    start = start+ap1length
    apertures = struct.unpack(apertureformat2, fileContent[start:start+ap2length])
    start = start + ap2length
    logging.debug(f"apertures (nr, radius): {apertures}")

    ############################### Stars ###############################################

    # read nr of stars
    nrstars = struct.unpack("i", fileContent[start:start+4])[0]
    logging.debug(f"Number of stars: {nrstars}")
    start = start+4
    # read stars
    starformat = "<2i5d"
    starsize = struct.calcsize(starformat)
    stars = []
    Star = namedtuple('Star', 'id, ref_id, x, y, skymed, skysig, fwhm')
    for _ in range(nrstars):
        stars.append(Star._make(struct.unpack(starformat, fileContent[start:start+starsize])))
        start = start + starsize
    logging.debug(f"Read {len(stars)} stars, first 10: {stars[0:10]}")

    ############################## Data
    dataformat = "<iii"
    dataformatsize = struct.calcsize(dataformat)
    starsleft = (len(fileContent) - start)/(dataformatsize*napertures)
    logging.debug(f"starting at {start} so we can read {starsleft}")
    StarData = namedtuple('StarData', 'mag, err, code')
    INT_MAX = 2147483647
    stardata = []
    for _ in range(int(starsleft)):
        result = []
        for apidx in range(napertures):
            if only_apertureidx != -1 and apidx != only_apertureidx:
                result.append(None)
            else:
                mag, err, code = struct.unpack(dataformat, fileContent[start:start+dataformatsize])
                if mag != INT_MAX:
                    mag = mag / 0x1000000
                else:
                    mag = sys.float_info.max
                if err != INT_MAX:
                    err = err / 0x1000000
                if code != 0 and mag != sys.float_info.max:
                    logging.debug(f"Error code: {code} {mag}") # see cmpack_common.h, enum CmpackError
                result.append(StarData._make((mag, err, code)))
            start = start+dataformatsize
        stardata.append(result)

    logging.debug(f"Read {len(stardata)} data for stars, first star: {stardata[0]}")
    logging.debug(len(stardata[0]))
    logging.debug(f"size in memory is: {(sys.getsizeof(stardata)+sys.getsizeof(stars))/1024}")
    return (photheader, apertures, nrstars, stars, stardata)

# percentage is between 0 and 1
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
    logging.debug(files)
    # pre process
    apertures = []
    with open(files[0], mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        photheader, apertures, nrstars , _, _ = read_pht_file(files[0], fileContent)
        apertures = convert_to_aperture_only(apertures)
    file.close()
    logging.info(f"Apertures: {apertures}, nr of files: {nrfiles}")
    collect = np.empty([len(apertures), nrstars, nrfiles, 2],dtype=float)
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

    print("Nr of stars for aperture 2:", len(collect[2]))
    print("Mb used:", collect.nbytes/1024/1024)
    print("Entry for First aperture, first star:", collect[0][0])
    print("Shape for collect:", collect.shape)
    stddevs = np.empty([len(apertures), nrstars])
    for apidx in range(len(collect)):
        for staridx in range(len(collect[apidx])):
            # print(apidx, staridx)
            stddevs[apidx, staridx] = collect[apidx][staridx].std()
    print(stddevs[0,0:20])
    print("Shape for stddevs:", stddevs.shape)
    print(np.argmin(stddevs[0]))
    # for idx in range(len(stddevs)):
    #     print(stddevs[idx].min(), stddevs[idx].max())
    #     print(idx, stddevs[idx].sum())
    median_erroradded = np.median(np.add(np.take(fwhm, 1, axis=1), np.take(fwhm, 2, axis=1)))
    median_multiply = np.median(np.take(fwhm, 1, axis=1))*1.75

    apertureidx = np.abs(apertures - median_multiply).argmin()
    print("FWHM median:", median_multiply, "aperture chosen is:", apertures[apertureidx])

    compstars_0 = np.argpartition(stddevs[apertureidx], 3)[:3]

    #compstar_0 = np.argmin(stddevs[apertureidx], axis=0)
    print("Compstars_0 with minimum stdev in the chosen aperture:", compstars_0)
    return stddevs, collect, apertures, apertureidx, fwhm, jd, (compstars_0+1).tolist()

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    read_pht_file(sys.argv[1], open(sys.argv[1],'rb').read())

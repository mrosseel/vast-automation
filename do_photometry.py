import os
import struct
from collections import namedtuple
import glob
import init
import sys
import numpy as np
import logging


def process_file(the_file, fileContent):
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
    napertures = struct.unpack(apertureformat1, fileContent[start:start+ap1length])[0]
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
    for star in range(nrstars):
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
    for stentry in range(int(starsleft)):
        result = []
        for apentry in range(napertures):
    #         try:
            mag, err, code = struct.unpack(dataformat, fileContent[start:start+dataformatsize])
            if mag != INT_MAX:
                mag = mag / 0x1000000
            else:
                mag = sys.float_info.max
            if err != INT_MAX:
                err = err / 0x1000000
            if code != 0 and mag != sys.float_info.max:
                logging.debug("Error code:", code, mag) # see cmpack_common.h, enum CmpackError
            result.append(StarData._make((mag, err, code)))
            start = start+dataformatsize
    #         except:
    #             logging.debug(stentry, apentry, start)
    #             logging.debug(dataformatsize)
        stardata.append(result)

    logging.debug(f"Read {len(stardata)} data for stars, first star: {stardata[0]}")
    logging.debug(len(stardata[0]))
    logging.debug(f"size in memory is: {(sys.getsizeof(stardata)+sys.getsizeof(stars))/1024}")
    return (photheader, apertures, nrstars, stars, stardata)

# percentage is between 0 and 1
def select_files(match_pattern, percentage=1):
    matched_files = glob.glob(init.matchedphotometrydir+match_pattern)
    desired_length = int(len(matched_files) * percentage)
    np.random.seed(42) # for the same percentage, we always get the same selection
    selected_files = np.random.choice(matched_files, size=desired_length, replace=False).tolist()
    return selected_files

def main():
    files = select_files('match*.pht', 0.05)
    nrfiles = len(files)
    logging.debug(files)
    # pre process
    apertures = []
    with open(files[0], mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        apertures = process_file(files[0], fileContent)[1][1::2]
    logging.info(f"Apertures: {apertures}")

    for fileidx, entry in enumerate(files):
        fileContent = None
        with open(entry, mode='rb') as file: # b is important -> binary
            fileContent = file.read()
            photheader, apertures, nrstars, stars, stardata = process_file(entry, fileContent)
            # logging.debug(f"the result is {result[0]}")
            collect = np.empty([len(apertures), nrstars, nrfiles],dtype=float)
            for staridx, starentry in enumerate(stars):
                for apidx, aperturedata in enumerate(stardata[staridx]):
                    collect[apidx][starentry.ref_id-1][fileidx] = aperturedata.mag
                # if index == 0:
                #     print(index, entry, stardata[index])
            #     collect[entry.id-1][0] = entry.id
            #     collect[entry.id-1][1] = entry.ref_id
            #     collect[entry.id-1][2] = entry.fwhm
    print(collect[0][0])



if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    main()

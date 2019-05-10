import struct
from collections import namedtuple
import sys
import numpy as np
import logging


# header
headerformat = "<28sxxxx4id70sdd70sH5Bxddddddddddddiiii3d70sdd70sdddddddd"
headerlength = struct.calcsize(headerformat)
unpack_headerformat = struct.Struct(headerformat).unpack
PhotHeader = namedtuple('PhotHeader', 'magic, revision, headerlen, width, height, jd, filter, exptime, ccdtemp, origin, year, mon, mday, hour, min, sec, \
                range0, range1, gain, rnoise, fwhm_exp, fwhm_mean, fwhm_err, threshold, sharpness0, sharpness1, roundness0, roundness1, matched, match_rstars, \
                match_istars, match_mstars, match_clip, offset0, offset1, objectdesig, objectra, objectdec, locationdesignation, \
                locationlon, locationlat, trafoxx, trafoxy, trafox0, trafoyx, trafoyy, trafoy0 \
                ')
# WCS
wcsformat= "<i"
wcssize = struct.calcsize(wcsformat)
unpack_wcsformat = struct.Struct(wcsformat).unpack
# apertures
apertureformat1 = "<i"
ap1length = struct.calcsize(apertureformat1)
unpack_apertureformat1 = struct.Struct(apertureformat1).unpack
# star
nrstarsformat = "i"
nrstarssize = struct.calcsize(nrstarsformat)
unpack_nrstarsformat= struct.Struct(nrstarsformat).unpack
starformat = "<2i5d"
starsize = struct.calcsize(starformat)
unpack_star =struct.Struct(starformat).unpack
Star = namedtuple('Star', 'id, ref_id, x, y, skymed, skysig, fwhm')
# data
dataformat = "<iii"
dataformatsize = struct.calcsize(dataformat)
unpack_data =struct.Struct(dataformat).unpack
StarData = namedtuple('StarData', 'mag, err, code')
INT_MAX = 2147483647


np_dtype = np.dtype([('id', np.int), ('ref_id', np.int), ('x', np.float), ('y', np.float), ('skymed', np.float),
                     ('skysig', np.float), ('fwhm', np.float)])
# add check for magic string?
MAGIC_STRING = "C-Munipack photometry file\r\n" # length is 28

# if only_apertureidx is provided, the returned stardata is a flat arrray, otherwise it's 2d (stars, apertures)
def read_pht_file(the_file, fileContent, only_apertureidx: int =-1, read_stars=True, check=True):
    assert type(only_apertureidx) == int
    logging.debug(f"The file:{the_file}")
    ############################### HEADER ###############################################

    logging.debug(f"header length is {headerlength}")
    unp = unpack_headerformat(fileContent[:headerlength])
    photheader = PhotHeader._make(unp)
    if check:
        assert bytes(MAGIC_STRING, encoding='utf-8') == photheader.magic
    for name, entry in photheader._asdict().items():
        logging.debug(f"name entries: {name}: {entry}")
    ############################### WCS ###############################################
    wcsdatalength = unpack_wcsformat(fileContent[headerlength:headerlength+wcssize])[0]
    wcslength = headerlength + wcssize + wcsdatalength
    logging.debug(f"wcs data length: {wcsdatalength}")
    ############################### APERTURES ###############################################
    # typedef struct _CmpackPhtAperture
    # {
    # 	int id;								/**< Aperture identifier */
    # 	double radius;						/**< Radius in pixels */
    # } CmpackPhtAperture;

    start = wcslength
    napertures = unpack_apertureformat1(fileContent[start:start+ap1length])[0]
    logging.debug(f"napertures: {napertures}")
    apertureformat2 = "<" + napertures * "id"
    ap2length = struct.calcsize(apertureformat2)
    start += ap1length
    apertures = struct.unpack(apertureformat2, fileContent[start:start+ap2length])
    start += ap2length
    logging.debug(f"apertures (nr, radius): {apertures}")

    ############################### Stars ###############################################

    # read nr of stars
    nrstars = unpack_nrstarsformat(fileContent[start:start+4])[0]
    start = start+4
    logging.debug(f"Number of stars: {nrstars}")
    # read stars

    stars = []
    if read_stars:
        for _ in range(nrstars):
            stars.append(Star._make(unpack_star(fileContent[start:start+starsize])))
            start += starsize
        logging.debug(f"Read {len(stars)} stars, first 10: {stars[0:10]}")
    else:
        start += starsize*nrstars

    ############################## Data
    # starsleft = (len(fileContent) - start)/(dataformatsize*napertures)
    logging.debug(f"starting at {start} so we can read {nrstars}")
    stardata = []
    if only_apertureidx != -1:
        start += only_apertureidx*dataformatsize

    for idx in range(nrstars):
        result = []
        if only_apertureidx != -1:
            mag, err, code = unpack_data(fileContent[start:start+dataformatsize])
            result = StarData._make((getConvertedStarData(mag, err, code, idx, the_file)))
            start += napertures*dataformatsize
        else:
            # logging.debug(f"raw bytes: {fileContent[start:start+dataformatsize]}")
            for apidx in range(napertures):
                mag, err, code = unpack_data(fileContent[start:start+dataformatsize])
                result.append(StarData._make((getConvertedStarData(mag, err, code, idx, the_file))))
                start += dataformatsize
        stardata.append(result)

    logging.debug(f"Read {len(stardata)} data for stars, first star: {stardata[0]}")
    logging.debug(len(stardata[0]))
    logging.debug(f"size in memory is: {(sys.getsizeof(stardata)+sys.getsizeof(stars))/1024}")
    return (photheader, apertures, nrstars, stars, stardata)

# take magnitude and error, and convert them to NAN if it's nonsense
def getConvertedStarData(mag, err, code, idx, the_file):
    #print(f"Converting: {mag}, {err}, {code}")
    if mag != INT_MAX:
        mag = mag / 0x1000000
    else:
        mag = np.nan
    if err != INT_MAX:
        err = err / 0x1000000
    else:
        err = np.nan
    if code != 0:
        logging.debug(f"Error code: {code} mag:{mag}, err:{err}, star_0:{idx}, file:{the_file}") # see cmpack_common.h, enum CmpackError
    #print(f"Converted: {mag}, {err}, {code}")
    return mag, err, code

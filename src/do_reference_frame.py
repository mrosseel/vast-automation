# https://github.com/deprecated/fits2image
import sys
import numpy as np
try:
    from astropy.io import fits
except ImportError:
    import pyfits as fits
from PIL import Image
from io import BytesIO
import tqdm
import multiprocessing as mp
from functools import partial
import glob
import init

def fits_to_size(fitsfilename, vmin, vmax):
    return fitsdata_to_size(get_fits_data(fitsfilename, vmin, vmax), fitsfilename)

def fits_to_jpeg(fitsfilename, vmin, vmax):
    return fitsdata_to_jpeg(get_fits_data(fitsfilename, vmin, vmax), fitsfilename)

def get_fits_data(fitsfilename, vmin, vmax):
    # Try to read data from first HDU in fits file
    data = fits.open(fitsfilename)[0].data
    # If nothing is there try the second one
    if data is None:
        data = fits.open(fitsfilename)[1].data

    maxdata = 0
    mindata = sys.maxsize
    #print(np.min(data), np.max(data))

    # Clip data to brightness limits
    data[data > vmax] = vmax
    data[data < vmin] = vmin
    # Scale data to range [0, 1]
    data = (data - vmin)/(vmax - vmin)
    # Convert to 8-bit integer
    data = (255*data).astype(np.uint8)
    # Invert y axis
    data = data[::-1, :]
    return data


def fitsdata_to_size(data, fitsfilename):
    # Create image from data array and save as jpg
    image = Image.fromarray(data, 'L')
    img_file = BytesIO()
    image.save(img_file, 'jpeg')
    image_file_size = img_file.tell()
    return (fitsfilename, image_file_size)

def fitsdata_to_jpeg(data, fitsfilename):
    # Create image from data array and save as jpg
    image = Image.fromarray(data, 'L')
    imagename = fitsfilename.replace('.fts', '.jpg')
    image.save(imagename)
    return imagename

# limit is not used at the moment
def runit(fitsdir, limit):
    # between 0 and 65535
    _vmin=50
    _vmax=5000
    result = []
    fitslist = glob.glob(fitsdir+'/*.fts')
    pool = mp.Pool(init.nr_threads, maxtasksperchild=100)
    func = partial(fits_to_size, vmin=_vmin, vmax=_vmax)
    print("Calculating fits quality for", len(fitslist), " images in", fitsdir)
    for resulttuple in tqdm.tqdm(pool.imap_unordered(func, fitslist, 10), total=len(fitslist)):
        result.append(resulttuple)
    bestfile = ''
    bestsize = 0
    for name, size in result:
        if size > bestsize:
            bestsize = size
            bestfile = name
    jpegfile = fits_to_jpeg(bestfile, _vmin, _vmax)
    print("result:", bestfile, bestsize, jpegfile)
    return (bestfile, bestsize, jpegfile)

if __name__ == '__main__':
    runit(init.convfitsdir)

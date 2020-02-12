import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import re
from pathlib import Path
import matplotlib.animation as animation
from tqdm import tqdm
from scipy import ndimage
from matplotlib.colors import LogNorm
import reading
from collections import namedtuple
from datetime import datetime
import logging
import argparse
import os
import utils

ImageRecord = namedtuple('ImageRecord', 'jd, x, y, file, rotation')
im = []
fitsdir = ""
crop = 0
fig = plt.figure(figsize=(15, 15), dpi=80, facecolor='w', edgecolor='k')
ax = plt.axes()
pbar = None
padding = 200
fullcrop = 0


def animate(vastdir: str, resultdir: str, afitsdir: str, starid: int, acrop: int):
    global fitsdir
    global crop
    global im
    global pbar
    global fullcrop
    crop = acrop
    fitsdir = afitsdir
    im = ax.imshow(np.zeros(crop * crop).reshape(crop, crop), cmap='gray', origin='lower', vmin=0, vmax=2500)
    border = 20
    fullcrop = crop + border
    dpi = 100

    # read the lightcurve of the chosen star
    lightcurvefile = f'out{starid:05}.dat'
    logging.info(f"starid {starid}, file is {lightcurvefile}")
    df = reading.read_lightcurve_vast(Path(vastdir, lightcurvefile))
    image_records = []
    rotation_dict = extract_frame_rotation_dict(vastdir)
    for index, row in df.iterrows():
        filename = Path(row['file']).name
        image_records.append(ImageRecord(float(row['JD']), round(row['X']), round(row['Y']), filename, float(rotation_dict[filename])))

    logging.info(f"rotation dict has {len(rotation_dict)} entries")
    logging.info(f"imagerecords has {len(image_records)} entries")
    sorted_images = sorted(image_records, key=lambda x: x.jd)

    pbar = tqdm(total=len(sorted_images))
    ani = animation.FuncAnimation(fig, update_img, frames=sorted_images, interval=30, init_func=init)
    writer = animation.writers['ffmpeg'](fps=30)

    ani.save(Path(resultdir, f'demo.mp4'), writer=writer, dpi=dpi)
    pbar.close()


# get the mapping fits file -> rotation
def extract_frame_rotation_dict(vastdir) -> float:
    filename = Path(vastdir, 'vast_image_details.log')
    the_regex = re.compile(r'^.*rotation=\s*([0-9,.,-]+).*\s+(.+)$')
    rotation_dict = {}
    with open(filename, 'r') as infile:
        for line in infile:
            thesearch = the_regex.search(line)
            if thesearch:
                path = Path(thesearch.group(2))
                rotation_dict[path.name] = float(thesearch.group(1).strip())
    return rotation_dict


# make a dict with mapping fits file -> image00007.cat
def extract_image_catalog(vastdir) -> float:
    filename = Path(vastdir, 'vast_images_catalogs.log')
    the_regex = re.compile(r'^(.*) (.*)$')
    catalog_dict = {}
    with open(filename, 'r') as infile:
        for line in infile:
            thesearch = the_regex.search(line)
            if thesearch:
                path = Path(thesearch.group(2))
                catalog_dict[path.name] = thesearch.group(1)
    return catalog_dict


def crop_center(img, cropx, cropy):
    y, x = img.shape
    startx = x // 2 - (cropx // 2)
    starty = y // 2 - (cropy // 2)
    return img[starty:starty + cropy, startx:startx + cropx]


def init():
    pass


def update_img(record: ImageRecord):
    hdulist = fits.open(Path(fitsdir, record.file))
    data = hdulist[0].data.astype(float)
    shapex, shapey = hdulist[0].shape
    backgr = data.mean()
    data = data.reshape(shapex, shapey)
    data = np.pad(data, (padding, padding), 'constant', constant_values=(backgr, backgr))
    cropdata = data[record.y - fullcrop + padding:record.y + fullcrop + padding,
               record.x - fullcrop + padding:record.x + fullcrop + padding]
    rotcrop = ndimage.interpolation.rotate(cropdata, record.rotation)
    rotcrop = crop_center(rotcrop, crop, crop)
    #fig = plt.figure(figsize=(15, 15), dpi=80, facecolor='w', edgecolor='k')

    # rotx, roty = rotcrop.shape
    # target_app = CircularAperture((rotx // 2, roty // 2), r=5.)
    # target_app.plot()
    # data[record.y-cropsize:record.y+cropsize,record.x-cropsize:record.x+cropsize] = 0
    # ax.imshow(rotcrop, cmap='gray', origin='lower', norm=LogNorm())
    pbar.update(1)
    im.set_data(rotcrop)
    ax.set_title(f"JD: {record.jd}")
    # plt.suptitle(f"Nr {idx}, rotation is {record.rotation}, shape is {rotcrop.shape}", size=16)


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description='munipack automation cli')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (usually the vast dir)",
                        nargs='?', required=True)
    parser.add_argument('-r', '--resultdir',
                        help="The directory where all results will be written",
                        nargs='?', required=True)
    parser.add_argument('--fitsdir', help="The dir where the fits are, only needed to plate-solve the reference frame",
                        required=False)
    parser.add_argument('-s', '--star',
                        help="The star from which we want to make a movie",
                        nargs='?', required=True)
    parser.add_argument('-c', '--crop',
                        help="The width of",
                        nargs='?', required=True)
    # parser.add_argument('-s', '--start',
    #                     help="The start Julian Date",
    #                     nargs='?', required=True)
    # parser.add_argument('-e', '--end',
    #                     help="The end Julian Date",
    #                     nargs='?', required=True)
    parser.add_argument('-x', '--verbose', help="Set logging to debug mode", action="store_true")
    args = parser.parse_args()
    datadir = utils.add_trailing_slash(args.datadir)
    datenow = datetime.now()
    resultdir = Path(args.resultdir)
    filehandler = f"{datadir}vastlog-{datenow:%Y%M%d-%H_%M_%S}.log"
    fh = logging.FileHandler(filehandler)
    fh.setLevel(logging.INFO)
    # add the handlers to the logger
    logger.addHandler(fh)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        fh.setLevel(logging.DEBUG)

    assert os.path.exists(args.datadir), "datadir does not exist"
    # assert os.path.exists(args.resultdir), "resultdir does not exist" ==> this dir is created
    assert os.path.exists(args.resultdir), "resultdir does not exist"
    assert os.path.exists(args.fitsdir), "fitsdir does not exist"
    animate(datadir, args.resultdir, args.fitsdir, int(args.star), int(args.crop))

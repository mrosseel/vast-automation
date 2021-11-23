from typing import List, Dict, Tuple
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from functools import partial
from photutils import aperture_photometry, CircularAperture
from astropy.coordinates import match_coordinates_sky
from pathlib import Path
import main_vast
import utils
import os
from reading import ImageRecord
import do_calibration
import utils_sd
import subprocess
import logging
import time
import reading
import argparse
from datetime import datetime
import random
import toml
from collections import namedtuple

from star_description import StarDescription
from ucac4 import UCAC4, MinimalStarTuple

padding = 0
dpi = 600
NEIGHBOUR_TEXT_SIZE = 4
UCAC4_TEXT_SIZE = 3
RefFrame = namedtuple("RefFrame", "ref_jd path_to_solved path_to_reference_frame")
ucac4 = UCAC4()


def inspect(vastdir, resultdir, fitsdir, apikey, stars):
    """generate image+txt file to inspect a certain star on wrong ucac/localid nrs """
    reading.trash_and_recreate_dir(Path(resultdir) / "inspect")
    ref_jd, _, _, reference_frame = reading.extract_reference_frame(vastdir)
    _, shapex, shapey = reading.get_fits_data(Path(fitsdir, Path(reference_frame).name))
    refframes: List[RefFrame] = [
        RefFrame(
            ref_jd,
            Path(vastdir) / "new-image.fits",
            Path(fitsdir, Path(reference_frame).name),
        )
    ]

    # construct star descriptions
    sds = utils_sd.construct_star_descriptions(vastdir, None)
    star_catalog = do_calibration.create_star_descriptions_catalog(sds)

    if args.radecs:
        main_vast.read_and_tag_radec(args.radecs, sds)
    if args.localids:
        stardict = main_vast.get_localid_to_sd_dict(sds)
        main_vast.read_and_tag_localid(args.localids, stardict)
    if args.stars:
        stars = sorted(list(map(lambda x: int(x), stars)))
    else:
        stars = sorted([x.local_id for x in sds if x.get_metadata("SELECTEDTAG")])
    for starid in stars:
        process(
            vastdir,
            resultdir,
            fitsdir,
            apikey,
            shapex,
            shapey,
            starid,
            ref_jd,
            reference_frame,
            sds,
            star_catalog,
            refframes,
        )


def process(
    vastdir,
    resultdir,
    fitsdir,
    apikey,
    shapex,
    shapey,
    starid,
    ref_jd,
    reference_frame,
    sds,
    star_catalog,
    refframes: List[RefFrame],
):
    logging.info(f"Processing star {starid}")
    # getting the image records and rotation dict for star X
    image_records, rotation_dict = reading.get_star_jd_xy_rot(starid, vastdir)
    logging.debug(f"rotation dict has {len(rotation_dict)} entries")
    logging.debug(f"imagerecords has {len(image_records)} entries")
    df = pd.DataFrame(image_records)
    jddict = reading.jd_to_fitsfile_dict(vastdir)
    platesolved_file = None

    # search a refframe where the star is present
    for frame in refframes:
        # check if reference frame is good enough
        refdf = df[df.jd == float(ref_jd)]
        if (
            len(refdf) == 1
            and (abs(shapex / 2 - refdf.iloc[0].x) < shapex / 2 * 0.95)
            and (abs(shapey / 2 - refdf.iloc[0].y) < shapey / 2 * 0.95)
        ):
            logging.info(
                f"Choosing refframe {frame.path_to_reference_frame.name} because star {starid} is on it."
            )
            refrow = refdf.iloc[0]
            chosen_record = ImageRecord(refrow.jd, refrow.x, refrow.y, refrow.file, 0.0)
            chosen_fits_fullpath = frame.path_to_reference_frame
            platesolved_file = frame.path_to_solved
            logging.info("Found platesolved file.")
            break
    # couldn't find a match in the existing refframes, let's find a new one
    if not platesolved_file:
        logging.info("No existing frames can be used for this star.")
        df = df[(df.x > 0) & (df.x < shapex) & (df.y > 0) & (df.y < shapey)]
        df["borderdistance"] = abs(shapex / 2 - df.x) + abs(shapey / 2 - df.y)
        logging.debug(df["borderdistance"].idxmin())
        logging.debug(df.iloc[df["borderdistance"].idxmin()])
        chosen_row = df.iloc[df["borderdistance"].idxmin()]
        chosen_jd = chosen_row.jd
        chosen_record = ImageRecord(
            chosen_jd, chosen_row.x, chosen_row.y, chosen_row.file, chosen_row.rotation
        )
        chosen_fits_fullpath = Path(fitsdir, jddict[chosen_jd])
        platesolved_file = Path(vastdir) / f"platesolve_{starid}.fits"
        if os.path.isfile(platesolved_file):
            logging.info(f"Found old platesolved file on disk: {platesolved_file}")
        else:
            logging.info(f"Platesolving {platesolved_file}")
            plate_solve(apikey, chosen_fits_fullpath, platesolved_file)
        refframes.append(RefFrame(chosen_jd, platesolved_file, chosen_fits_fullpath))
    chosen_rotation = rotation_dict[chosen_fits_fullpath.name]
    logging.info(
        f"platesolved file: {platesolved_file}, chosen fits: {chosen_fits_fullpath}, "
        f"rotation: {chosen_rotation}"
    )
    get_closest_stars(sds, star_catalog, chosen_record, resultdir, platesolved_file, starid)

# Get the 10 closest neighbour stars
def get_closest_stars(sds, star_catalog, chosen_record, resultdir, platesolved_file, starid):
    sd_dict = utils.get_localid_to_sd_dict(sds)
    chosen_star_sd = sd_dict[starid]
    neighbours = []
    for neigh in range(2, 22):
        idx, d2d, _ = match_coordinates_sky(
            chosen_star_sd.coords, star_catalog, nthneighbor=neigh
        )
        neighbours.append(sds[idx])
    ucac4.add_sd_metadatas(neighbours)
    ucac4.add_sd_metadatas([chosen_star_sd])
    update_img(chosen_star_sd, chosen_record, neighbours, resultdir, platesolved_file)


# platesolve the best frame for this star
def plate_solve(apikey, chosen_fits_fullpath, output_file):
    subprocess.Popen(
        f"python3 ./src/astrometry_api.py --apikey={apikey} "
        f"--upload={chosen_fits_fullpath} --newfits={output_file} --private --no_commercial",
        shell=True,
    )
    while not os.path.isfile(output_file):
        logging.info(f"Waiting for the astrometry.net plate solve...")
        time.sleep(10)


#
def update_img(
    star: StarDescription,
    record: ImageRecord,
    neighbours: List[StarDescription],
    resultdir: str,
    platesolved_file: str,
):
    resultlines = []
    fig = plt.figure(figsize=(36, 32), dpi=dpi, facecolor="w", edgecolor="k")
    wcs = do_calibration.get_wcs(platesolved_file)
    data, shapex, shapey = reading.get_fits_data(platesolved_file)
    backgr = data.mean()
    data = data.reshape(shapex, shapey)
    data = np.pad(
        data, (padding, padding), "constant", constant_values=(backgr, backgr)
    )
    starxy = SkyCoord.to_pixel(star.coords, wcs=wcs, origin=0)

    # add main target
    add_circle(record.x, record.y, 3, "b")
    startoml = load_toml(star, resultdir)
    star.vmag = startoml["vmag"]
    resultlines.append(log_star(star, -1))
    random_offset = False
    offset1 = 70
    offset2 = 10

    # add neighbours
    for idx, nstar in enumerate(neighbours):
        add_pixels(nstar, wcs, 0)
        add_circle(nstar.xpos, nstar.ypos, 4, "g")
        if random_offset:
            xrandoffset = random.randint(3, 4)
            yrandoffset = random.randint(2, 3)
            xsignrand = random.choice([-1.0, 1.0])
            ysignrand = random.choice([-1.0, 1.0])
            offset1 = xsignrand * xrandoffset
            offset2 = ysignrand * yrandoffset
        # https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.arrow.html
        plt.annotate(
            f"{idx}",
            xy=(round(nstar.xpos), round(nstar.ypos)),
            xycoords="data",
            xytext=(offset1, offset2),
            textcoords="offset points",
            size=NEIGHBOUR_TEXT_SIZE,
            color="red",
            arrowprops=dict(arrowstyle="-", color="grey", alpha=0.2),
        )
        resultlines.append(log_star(nstar, idx))

    # loading and painting ucac stars
    radius = 0.08
    ucac_stars: List[MinimalStarTuple] = ucac4.get_region_minimal_star_tuples(
        star.coords.ra.deg, star.coords.dec.deg, radius
    )
    logging.info(f"Looping on {len(ucac_stars)} UCAC4 stars")
    for ucac_star in ucac_stars:
        coord = SkyCoord(ucac_star.ra, ucac_star.dec, unit="deg")
        if star.coords.separation(coord).degree > radius:
            continue
        xy = SkyCoord.to_pixel(coord, wcs=wcs, origin=0)
        x, y = round(xy[0].item(0)), round(xy[1].item(0))
        add_circle(x, y, 2, "c")
        plt.annotate(
            f"{ucac_star.id[-3:]}",
            xy=(x, y),
            xycoords="data",
            xytext=(2, 2),
            textcoords="offset points",
            size=UCAC4_TEXT_SIZE,
            arrowprops=dict(arrowstyle="-", color="grey", alpha=0.2),
        )

    median = np.median(data)
    #     data = ndimage.interpolation.rotate(data, record.rotation)

    plt.imshow(data, cmap="gray_r", origin="lower", vmin=0, vmax=min(median * 5, 65536))
    starui = utils.get_star_or_catalog_name(star)
    save_inspect_image = Path(
        Path(resultdir) / "inspect",
        f"{starui.filename_no_ext}_{Path(resultdir).name}_inspect.png",
    )
    save_inspect_txt = Path(
        Path(resultdir) / "inspect",
        f"{starui.filename_no_ext}_{Path(resultdir).name}_inspect.txt",
    )
    fig.savefig(save_inspect_image)
    logging.info(f"Saved file as {save_inspect_image}.")
    write_file(star, save_inspect_txt, resultlines)
    plt.close(fig)
    plt.clf()


# UNUSED ATM
# returns cropped image and coordinate transformer
def crop_around(img, aroundx: float, aroundy: float, halfwidth: int):
    y, x = img.shape
    # startx = x // 2 - (cropx // 2)
    startx = max(0, int(aroundx) - halfwidth)
    # starty = y // 2 - (cropy // 2)
    starty = max(0, int(aroundy) - halfwidth)
    translater = partial(
        coord_translate, centerx=aroundx, centery=aroundy, halfwidth=halfwidth
    )
    return (
        img[starty : min(y, starty + halfwidth), startx : min(x, startx + halfwidth)],
        translater,
    )


# UNUSED ATM and not correct
def coord_translate(xin, yin, centerx, centery, halfwidth: int):
    return xin - centerx + halfwidth, yin - centery + halfwidth


def load_toml(star, resultdir):
    starui: utils.StarUI = utils.get_star_or_catalog_name(star)
    logging.info(f"starui: {starui}")
    txt_path = (
        Path(resultdir)
        / "phase_selected"
        / "txt"
        / f"{starui.filename_no_suff_no_ext}.txt"
    )
    if not os.path.isfile(txt_path):
        txt_path = (
            Path(resultdir) / "phase_selected" / "txt" / f"{star.local_id:05}.txt"
        )
    if not os.path.isfile(txt_path):
        txt_path = (
            Path(resultdir) / "phase_candidates" / "txt" / f"{star.local_id:05}.txt"
        )
    logging.debug("txtpath", txt_path)
    try:
        parsed_toml = toml.load(txt_path)
    except FileNotFoundError:
        logging.error(f"Could not load txt file with phase information from {txt_path}")
    return parsed_toml


def write_file(star: StarDescription, outputfile: str, lines: List[str]):
    logging.debug(f"Writing toml to {outputfile}")
    with open(outputfile, "w") as outfile:
        outfile.write("\n".join(lines))


def log_star(star, idx):
    star_ucac4 = star.get_metadata("UCAC4")
    starmag = f" mag: {star.vmag:.3f}" if star.vmag is not None else ""
    logline = (
        f"IDX: {idx}, Star: {star.local_id}{starmag} with ucac4: {star_ucac4.name}, "
        f"vmag: {star_ucac4.vmag}, coords: {star_ucac4.coords}"
    )
    logging.info(logline)
    return logline


def add_circle(xpos, ypos, radius: float, color: str):
    target_app = CircularAperture((xpos, ypos), r=radius)
    target_app.plot(color=color, alpha=0.7)


# kinda duplicate, move to util or something
def add_pixels(star, wcs, offset):
    star_coord = star.coords
    xy = SkyCoord.to_pixel(star_coord, wcs=wcs, origin=0, mode="all")
    x, y = round(xy[0].item(0)), round(xy[1].item(0))
    star.xpos, star.ypos = x, y


if __name__ == "__main__":
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description="munipack automation cli")
    parser.add_argument(
        "-d",
        "--datadir",
        help="The directory where the data can be found (usually the vast dir)",
        nargs="?",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--resultdir",
        help="The directory where all results will be written",
        nargs="?",
        required=True,
    )
    parser.add_argument(
        "--fitsdir",
        help="The dir where the fits are, only needed to plate-solve the reference frame",
        required=True,
    )
    parser.add_argument(
        "-a",
        "--radecs",
        help="Supply a file to identify stars known to you by RA/DEC",
        required=False,
    )
    parser.add_argument(
        "-l",
        "--localids",
        help="Load a file local ids, these ids will be used for field charts/reporting",
    )
    parser.add_argument(
        "--apikey",
        "-k",
        dest="apikey",
        help="API key for Astrometry.net web service; if not given will check AN_API_KEY environment variable",
    )
    parser.add_argument(
        "-s", "--stars", help="List the star id's to plot", nargs="+", required=False
    )
    parser.add_argument(
        "-x", "--verbose", help="Set logging to debug mode", action="store_true"
    )
    args = parser.parse_args()
    datadir = utils.add_trailing_slash(args.datadir)
    datenow = datetime.now()
    filehandler = f"{datadir}vastlog-{datenow:%Y%M%d-%H_%M_%S}.log"
    fh = logging.FileHandler(filehandler)
    fh.setLevel(logging.INFO)
    # add the handlers to the logger
    logger.addHandler(fh)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        fh.setLevel(logging.DEBUG)

    assert os.path.exists(args.datadir), "datadir does not exist"
    assert os.path.exists(args.resultdir), "resultdir does not exist"
    assert os.path.exists(args.fitsdir), "fitsdir does not exist"
    inspect(datadir, args.resultdir, args.fitsdir, args.apikey, args.stars)

from astropy.coordinates import SkyCoord
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils import aperture_photometry, CircularAperture
import numpy as np
import do_calibration
import logging
from reading import trash_and_recreate_dir
import argparse
from typing import List, Tuple
from star_description import StarDescription
import random
StarDescriptionList = List[StarDescription]

PADDING = 200
Shape = Tuple[int, int]

def set_local_id_label(star_descriptions):
    for star_descr in star_descriptions:
        star_descr.label = star_descr.local_id
    return star_descriptions

def set_aavso_id_label(star_descriptions):
    for star_descr in star_descriptions:
        star_descr.label = star_descr.aavso_id
    return star_descriptions

def set_custom_label(star_descriptions, label):
    for index, star_descr in enumerate(star_descriptions):
        star_descr.label = label if not isinstance(label, list) else label[index]
    return star_descriptions

def add_pixels(results, wcs, offset):
    for star in results:
        star_coord = star.coords
        xy = SkyCoord.to_pixel(star_coord, wcs=wcs, origin=0)
        x = xy[0].item(0)
        y = xy[1].item(0)
        star.xpos = x + offset
        star.ypos = y + offset
    return results


def mirror_offset_transform(pos: int, shape: Tuple, shapeidx: int, offset: int = 0):
    return shape[shapeidx] - (pos+offset)


def offset_transform(pos: int, shape: Tuple, shapeidx: int, offset: int = 0):
    return pos+offset


def plot_it(big_green: StarDescriptionList, small_red: StarDescriptionList, fits_file: str, wcs, title,
            padding: int = PADDING, plot_fits: bool = True, random_red=False,
            xpos_transform=offset_transform, ypos_transform=offset_transform):
    fig, data = get_plot_with_background(fits_file, padding, title, plot_fits)
    datashape = data.shape
    logging.info("plotting {} green and {} red circles.".format(len(big_green), len(small_red)))
    big_green = add_pixels(big_green, wcs, PADDING)
    small_red = add_pixels(small_red, wcs, PADDING)
    big_green_positions = [(o.xpos, o.ypos) for o in big_green]
    small_red_positions = [(o.xpos, o.ypos) for o in small_red]
    if len(big_green_positions) > 0:
        big_green_apps = CircularAperture(big_green_positions, r=10.)
        big_green_apps.plot(color='green', lw=1.5, alpha=0.5)
    if len(small_red_positions) > 0:
        small_red_apps = CircularAperture(small_red_positions, r=5.)
        small_red_apps.plot(color='red', lw=1.5, alpha=0.5)
    random.seed(42)

    # plot background fits image if one is provided
    # target_app.plot(color='blue', lw=1.5, alpha=0.5)
    #to_plot = results
    def annotate_it(star_descriptions, offset1, offset2, random_offset=False, size=16):
        for stardescr in star_descriptions:
            xpos = xpos_transform(stardescr.xpos, datashape, 0)
            ypos = ypos_transform(stardescr.ypos, datashape, 1)
            logging.debug(f"Plotting {stardescr.label} {xpos} {ypos}")
            if random_offset:
                randoffset = random.randint(10,20)
                xsignrand = random.choice([-1.0, 1.0])
                ysignrand = random.choice([-1.0, 1.0])
                offset1 = xsignrand*randoffset
                offset2 = ysignrand*randoffset
            plt.annotate('{}'.format(stardescr.label), xy=(round(xpos), round(ypos)), xycoords='data',
                         xytext=(offset1, offset2), textcoords='offset points', size=size,
                         arrowprops=dict(arrowstyle="->"))
    annotate_it(big_green, 0, -20, size=12)
    annotate_it(small_red, -10, 10, random_offset=random_red, size=10)
    return fig


#  plot_fits is false if no background needs to be plotted, in that case all zeros are used as data
def get_plot_with_background(fits_file: str, padding: int, title: str, plot_fits: bool = True):
    fig=plt.figure(figsize=(36, 32), dpi= 80, facecolor='w', edgecolor='k')
    plt.title(title, fontsize=40)
    hdulist = fits.open(fits_file)
    data = hdulist[0].data.astype(float)
    if not plot_fits:
        data = np.zeros(data.shape)
    data = np.pad(data, (padding, padding), 'constant', constant_values=(100, 100))
    plt.imshow(data, cmap='gray_r', origin='lower', vmin=0, vmax=2500)

    return fig, data


def save(fig, path):
    fig.savefig(path)
    plt.close(fig)


def run_standard_field_charts(star_descriptions: StarDescriptionList, wcs, fieldchartsdirs, reference_header):
    trash_and_recreate_dir(fieldchartsdirs)

    # setting the font size for titles/axes
    plt.rcParams.update({'axes.titlesize': 'large', 'axes.labelsize': 'large'})
    reference_fits_frame=reference_header
    SHOW_UPSILON = False

    # if SHOW_UPSILON:
    #     candidates = do_calibration.get_candidates(0.5)

    # TODO hand labeled stars
    # hand_candidates_descr = do_calibration.get_star_descriptions(init.wwcra_certain_candidates)
    all_stars_descr = star_descriptions

    # if SHOW_UPSILON:
    #     big_green = set_custom_label(comparison_  star_descr, 'comp')
    #     small_red = set_custom_label(apass_star_descr, [o.vmag for o in apass_star_descr])
    #     big_green = set_custom_label(vsx_star_descr, [o.match['catalog_dict']['name'] for o in vsx_star_descr])
    #     small_red = set_custom_label(hand_candidates_descr, [o.local_id for o in hand_candidates_descr])
    #     big_green = set_aavso_id_label(vsx_star_descr)
    #     small_red = set_local_id_label(hand_candidates_descr)

    # all stars get a blank label
    all_stars_labeled = set_custom_label(all_stars_descr, '')

    # vsx stars get their aavso id label
    vsx_descr = do_calibration.get_catalog_stars(star_descriptions, "VSX")
    vsx_labeled = set_aavso_id_label(vsx_descr)

    # candidate stars get their local id label
    candidate_descr = do_calibration.get_catalog_stars(star_descriptions, "CANDIDATE", exclude=["VSX"])
    candidate_labeled = set_local_id_label(candidate_descr)

    # starfile stars get their local id label
    starfile_descr = do_calibration.get_catalog_stars(star_descriptions, "STARFILE", exclude=["VSX"])
    starfile_labeled = set_local_id_label(starfile_descr)

    # default fields
    empty = []

    # field chart with all detections
    logging.info("Plotting field chart with all detected stars...")
    big_green = empty
    small_red = all_stars_labeled
    fig = plot_it(big_green, small_red, reference_fits_frame, wcs, "All detected stars", PADDING)
    save(fig, fieldchartsdirs + 'all_detections_{}_stars'.format(len(small_red)))

    # field chart with all vsx stars
    logging.info("Plotting field chart with all VSX variable stars...")
    big_green = vsx_labeled
    small_red = empty
    fig = plot_it(big_green, small_red, reference_fits_frame, wcs, "All VSX stars", PADDING)
    save(fig, fieldchartsdirs + 'vsx_stars_{}'.format(len(big_green)))

    # field chart with all vsx stars without the background
    logging.info("Plotting field chart with all VSX variable stars without reference field...")
    big_green = vsx_labeled
    small_red = empty
    fig = plot_it(big_green, small_red, reference_fits_frame, wcs, "VSX without background", PADDING, plot_fits=False)
    save(fig, fieldchartsdirs + 'vsx_stars_no_ref_{}'.format(len(big_green)))

    # field chart with only the background
    logging.info("Plotting field chart with only the reference field...")
    fig, _ = get_plot_with_background(reference_fits_frame, 0, "Reference frame")
    save(fig, fieldchartsdirs + 'only_ref')

    # field chart with all vsx stars + candidates
    logging.info("Plotting field chart with all VSX variable stars + candidate vars...")
    big_green = vsx_labeled
    small_red = candidate_labeled
    fig = plot_it(big_green, small_red, reference_fits_frame, wcs, "VSX stars + candidate stars", PADDING,
                  random_red=True)
    save(fig, fieldchartsdirs + 'vsx_{}_and_candidates_{}'.format(len(big_green), len(small_red)))

    # field chart with all vsx stars + starfile
    logging.info("Plotting field chart with all VSX variable stars + selected vars...")
    big_green = vsx_labeled
    small_red = starfile_labeled
    fig = plot_it(big_green, small_red, reference_fits_frame, wcs, "VSX stars + selected stars", PADDING,
                  random_red=True)
    save(fig, fieldchartsdirs + 'vsx_{}_and_selected_{}'.format(len(big_green), len(small_red)))

    # field charts for each individually selected starfile star
    big_green = vsx_labeled
    small_red = starfile_labeled
    for star in starfile_labeled:
        logging.info(f"Plotting field chart with all VSX variable stars + star {star.local_id}...")
        fig = plot_it(big_green, [star], reference_fits_frame, wcs, f"VSX stars + star {star.local_id}", PADDING,
                      random_red=False)
        save(fig, fieldchartsdirs + 'vsx_{}_and_selected_{}'.format(len(big_green), len(small_red)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='munipack automation field charts')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (fits in ./fits dir under the data dir",
                        nargs='?', required=True)
    parser.add_argument('-s', '--stars', help="List the star id's to plot", nargs='+')
    parser.add_argument('-n', '--novsx', help="Don't plot vsx stars", nargs='+')
    args = parser.parse_args()

from astropy.coordinates import SkyCoord
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils import aperture_photometry, CircularAperture
import numpy as np
from init_loader import init, settings
import do_calibration
import logging
from reading import trash_and_recreate_dir
import argparse
from typing import List
import utils
from star_description import StarDescription
StarDescriptionList = List[StarDescription]

PADDING = 200

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


def plot_it(big_green: StarDescriptionList, small_red: StarDescriptionList, fits_file, wcs, title, padding=PADDING):
    logging.info("plotting {} green and {} red circles.".format(len(big_green), len(small_red)))
    big_green = add_pixels(big_green, wcs, PADDING)
    small_red = add_pixels(small_red, wcs, PADDING)

    hdulist = fits.open(fits_file)
    data = hdulist[0].data.astype(float)
    data = np.pad(data, (padding, padding), 'constant', constant_values=(100, 100))
    fig=plt.figure(figsize=(36, 32), dpi= 80, facecolor='w', edgecolor='k')
    big_green_positions = ([o.xpos for o in big_green],[o.ypos for o in big_green])
    small_red_positions = ([o.xpos for o in small_red],[o.ypos for o in small_red])
    big_green_apps = CircularAperture(big_green_positions, r=10.)
    small_red_apps = CircularAperture(small_red_positions, r=5.)
    # target_app = CircularAperture(target_xy, r=20.)
    plt.title(title, fontsize=40)
    plt.imshow(data, cmap='gray_r', origin='lower', vmin=0, vmax=2500)
    big_green_apps.plot(color='green', lw=1.5, alpha=0.5)
    small_red_apps.plot(color='red', lw=1.5, alpha=0.5)
    # target_app.plot(color='blue', lw=1.5, alpha=0.5)
    #to_plot = results
    def annotate_it(star_descriptions, offset1, offset2, size=16):
        for stardescr in star_descriptions:
            plt.annotate('{}'.format(stardescr.label),
                         xy=(stardescr.xpos, stardescr.ypos), xycoords='data',
                         xytext=(offset1, offset2), textcoords='offset points', size=size, arrowprops=dict(arrowstyle="->"))
    annotate_it(big_green, -10, -20, size=10)
    annotate_it(small_red, -10, 10, size=12)
    return fig


def save(fig, path):
    fig.savefig(path)
    plt.close(fig)


def run_standard_field_charts(selected_star_descriptions: StarDescriptionList, wcs):
    trash_and_recreate_dir(settings.fieldchartsdirs)

    # setting the font size for titles/axes
    plt.rcParams.update({'axes.titlesize': 'large', 'axes.labelsize': 'large'})
    # reference_frame, reference_frame_index=reading.read_reference_frame()
    # reference_fits_frame=settings.convfitsdir+reference_frame
    reference_fits_frame=settings.reference_header
    SHOW_UPSILON = False

    # if SHOW_UPSILON:
    #     candidates = do_calibration.get_candidates(0.5)

    # TODO hand labeled stars
    # hand_candidates_descr = do_calibration.get_star_descriptions(init.wwcra_certain_candidates)
    all_stars_descr = do_calibration.get_star_descriptions()

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
    vsx_descr = [x for x in selected_star_descriptions if x.has_catalog('VSX')]
    vsx_labeled = set_aavso_id_label(vsx_descr)
    # other stars get their local id label
    hand_candidates_descr = [x for x in selected_star_descriptions if (x.has_catalog('SELECTED') and not x.has_catalog('VSX'))]
    hand_candidates_labeled = set_local_id_label(hand_candidates_descr)

    # default fields
    empty = []

    # field chart with all detections
    logging.info("Plotting field chart with all detected stars...")
    big_green = empty
    small_red = all_stars_labeled
    fig = plot_it(big_green, small_red, reference_fits_frame, wcs, "All detected stars", PADDING)
    save(fig, settings.fieldchartsdirs + 'all_detections_for_{}_stars'.format(len(small_red)))

    # field chart with all vsx stars
    logging.info("Plotting field chart with all VSX variable stars...")
    big_green = vsx_labeled
    small_red = empty
    fig = plot_it(big_green, small_red, reference_fits_frame, wcs, "All VSX variable stars", PADDING)
    save(fig, settings.fieldchartsdirs + 'all_vsx_stars_{}'.format(len(big_green)))

    # field chart with all vsx stars
    logging.info("Plotting field chart with all VSX variable stars + hand picked vars...")
    big_green = vsx_labeled
    small_red = hand_candidates_labeled
    fig = plot_it(big_green, small_red, reference_fits_frame, wcs, "VSX variable stars + hand picked stars", PADDING)
    save(fig, settings.fieldchartsdirs + 'all_vsx_stars_{}_hand_picked_{}'.format(len(big_green), len(small_red)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='munipack automation field charts')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (fits in ./fits dir under the data dir",
                        nargs='?', required=True)
    parser.add_argument('-s', '--stars', help="List the star id's to plot", nargs='+')
    parser.add_argument('-n', '--novsx', help="Don't plot vsx stars", nargs='+')
    args = parser.parse_args()

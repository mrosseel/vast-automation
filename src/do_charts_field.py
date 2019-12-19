from astropy.coordinates import SkyCoord
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils import CircularAperture
import numpy as np
import star_metadata
import utils
from comparison_stars import ComparisonStars
import logging
from reading import trash_and_recreate_dir
import argparse
from typing import List, Tuple
from star_description import StarDescription
import random
import do_compstars
import gc
import tqdm
StarDescriptionList = List[StarDescription]
gc.enable()
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


def set_custom_label(star_descriptions, label, strict=False):
    for index, star_descr in enumerate(star_descriptions):
        label_to_set = label if not isinstance(label, list) else label[index]
        if strict and (star_descr.label is not '' or None):
            logging.warning(f"Setting label to {label_to_set} but it has previous value {star_descr.label}")
        star_descr.label = label_to_set
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
    return shape[shapeidx] - (pos + offset)


def offset_transform(pos: int, shape: Tuple, shapeidx: int, offset: int = 0):
    return pos + offset


def annotate_it(star_descriptions, offset1, offset2, random_offset=False, size=16):
    for stardescr in star_descriptions:
        xpos = stardescr.xpos
        ypos = stardescr.ypos
        logging.debug(f"Field chart plotting label:{stardescr.label} x:{xpos} y:{ypos}")
        if random_offset:
            randoffset = random.randint(10, 20)
            xsignrand = random.choice([-1.0, 1.0])
            ysignrand = random.choice([-1.0, 1.0])
            offset1 = xsignrand * randoffset
            offset2 = ysignrand * randoffset
        plt.annotate(f'{stardescr.label}', xy=(round(xpos), round(ypos)), xycoords='data',
                     xytext=(offset1, offset2), textcoords='offset points', size=size,
                     arrowprops=dict(arrowstyle="->"))


'''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
RGB color; the keyword argument name must be a standard mpl colormap name.'''


def get_cmap(n, name='hsv'):
    return plt.cm.get_cmap(name, n)


def plot_it(star_lists: List[StarDescriptionList], sizes: List[float], random_offset: List[bool], fits_file: str, wcs, title,
            padding: int = PADDING, plot_fits: bool = True):
    fig, data = get_plot_with_background(fits_file, padding, title, plot_fits)
    logging.debug(f"plotting {[len(x) for x in star_lists]} stars per color")
    positions = []
    for stars in star_lists:
        stars = add_pixels(stars, wcs, PADDING)
        positions.append([(o.xpos, o.ypos) for o in stars])
    from itertools import cycle
    cycol = cycle('rgbcmk')

    for idx, pos in enumerate(positions):
        if len(pos) > 0:
            apps = CircularAperture(pos, r=sizes[idx])
            apps.plot(color=next(cycol), lw=1.5, alpha=0.5)
    random.seed(42)

    for idx, star in enumerate(star_lists):
        # annotate_it(big_green, 0, -20, size=12)
        annotate_it(star, -10, 15, random_offset=random_offset[idx], size=10)
    return fig


#  plot_fits is false if no background needs to be plotted, in that case all zeros are used as data
def get_plot_with_background(fits_file: str, padding: int, title: str, plot_fits: bool = True):
    fig = plt.figure(figsize=(36, 32), dpi=80, facecolor='w', edgecolor='k')
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


def run_standard_field_charts(star_descriptions: StarDescriptionList, wcs, fieldchartsdirs, reference_header,
                              comp_stars: ComparisonStars):
    trash_and_recreate_dir(fieldchartsdirs)

    # setting the font size for titles/axes
    plt.rcParams.update({'axes.titlesize': 'large', 'axes.labelsize': 'large'})
    reference_fits_frame = reference_header
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
    all_stars_no_label = set_custom_label(all_stars_descr, '')

    # vsx stars get their aavso id label
    vsx_descr = utils.get_stars_with_metadata(star_descriptions, "VSX")
    vsx_labeled = set_aavso_id_label(vsx_descr)

    # candidate stars get their local id label
    candidate_descr = utils.get_stars_with_metadata(star_descriptions, "CANDIDATE", exclude=["VSX"])
    candidate_labeled = set_local_id_label(candidate_descr)

    # starfile stars get their local id label
    starfile_descr = utils.get_stars_with_metadata(star_descriptions, "SITE")
    starfile_labeled = set_custom_label(starfile_descr, [x.get_metadata("SITE").our_name for x in starfile_descr])

    # owncatalog stars get their local id label
    owncatalog_descr = utils.get_stars_with_metadata(star_descriptions, "OWNCATALOG")
    owncatalog_labeled = set_custom_label(owncatalog_descr, [x.get_metadata("OWNCATALOG").name for x in owncatalog_descr])

    # field chart with all detections
    logging.info("Plotting field chart with all detected stars...")
    fig = plot_it([all_stars_no_label], [5.], [False], reference_fits_frame, wcs, "All detected stars", PADDING)
    save(fig, fieldchartsdirs + 'all_detections_{}_stars'.format(len(all_stars_no_label)))

    # field chart with all vsx stars
    logging.info("Plotting field chart with all VSX variable stars...")
    fig = plot_it([vsx_labeled], [10.], [False], reference_fits_frame, wcs, "All VSX stars", PADDING)
    save(fig, fieldchartsdirs + 'vsx_stars_{}'.format(len(vsx_labeled)))

    # field chart with all vsx stars without the background
    logging.info("Plotting field chart with all VSX variable stars without reference field...")
    fig = plot_it([vsx_labeled], [10.], [False], reference_fits_frame, wcs, "VSX without background", PADDING, plot_fits=False)
    save(fig, fieldchartsdirs + 'vsx_stars_no_ref_{}'.format(len(vsx_labeled)))

    # field chart with only the background
    logging.info("Plotting field chart with only the reference field...")
    fig, _ = get_plot_with_background(reference_fits_frame, 0, "Reference frame")
    save(fig, fieldchartsdirs + 'only_ref')

    # field chart with all vsx stars + candidates + owncatalog
    logging.info("Plotting field chart with all VSX variable stars + candidate vars...")
    fig = plot_it([vsx_labeled, candidate_labeled, owncatalog_labeled], [10., 5., 4.], [False, True, True],
                  reference_fits_frame, wcs, "VSX stars + candidate stars + own catalog", PADDING)
    save(fig, fieldchartsdirs + f'vsx_{len(vsx_labeled)}_and_candidates_{len(candidate_labeled)}')

    # field chart with all vsx stars + starfile + owncatalog
    logging.info(f"Plotting field chart with all VSX variable stars + {len(starfile_labeled)} selected vars...")
    fig = plot_it([vsx_labeled, starfile_labeled, owncatalog_labeled], [10., 5., 4.], [False, True, True],
                  reference_fits_frame, wcs, "VSX stars + selected stars + own catalog", PADDING)
    save(fig, fieldchartsdirs + 'vsx_{}_and_selected_{}'.format(len(vsx_labeled), len(starfile_labeled)))

    # compstars get their vmag as label
    # comp_stars_descr = comp_stars.star_descriptions
    # comp_stars_labeled = set_custom_label(comp_stars_descr, [x.vmag for x in comp_stars_descr])

    logging.info(f"Plotting field chart for each of the {len(starfile_labeled)} selected stars")
    # field charts for each individually selected starfile star
    for star in tqdm.tqdm(starfile_labeled):
        filtered_compstars_sds = do_compstars.get_star_compstars_from_catalog(star, comp_stars).star_descriptions
        compstars_labeled = set_custom_label(filtered_compstars_sds, [x.vmag for x in filtered_compstars_sds])
        filtered_compstars_sds = None

        fig = plot_it([[star], vsx_labeled, compstars_labeled], [10., 5., 3.], [True, False, True],
                      reference_fits_frame, wcs, f"VSX stars + star {star.local_id}", PADDING)
        save(fig, f"{fieldchartsdirs}vsx_and_star_{star.local_id:05}")
        gc.collect()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='munipack automation field charts')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (fits in ./fits dir under the data dir",
                        nargs='?', required=True)
    parser.add_argument('-s', '--stars', help="List the star id's to plot", nargs='+')
    parser.add_argument('-n', '--novsx', help="Don't plot vsx stars", nargs='+')
    args = parser.parse_args()

from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from photutils import CircularAperture
import numpy as np

import reading
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
        if strict and (star_descr.label != "" or None):
            logging.warning(
                f"Setting label to {label_to_set} but it has previous value {star_descr.label}"
            )
        star_descr.label = label_to_set
    return star_descriptions


def add_pixels(results, wcs, offset):
    for star in results:
        star_coord = star.coords
        xy = SkyCoord.to_pixel(star_coord, wcs=wcs, origin=0)
        x, y = round(xy[0].item(0)), round(xy[1].item(0))
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
        plt.annotate(
            f"{stardescr.label}",
            xy=(round(xpos), round(ypos)),
            xycoords="data",
            xytext=(offset1, offset2),
            textcoords="offset points",
            size=size,
            arrowprops=dict(arrowstyle="->"),
        )


"""Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
RGB color; the keyword argument name must be a standard mpl colormap name."""


def get_cmap(n, name="hsv"):
    return plt.cm.get_cmap(name, n)


def plot_it(
    star_lists: List[StarDescriptionList],
    sizes: List[float],
    random_offset: List[bool],
    fits_data: List[float],
    wcs,
    title,
    padding: int = PADDING,
    annotate=True,
):
    fig, data = get_plot_with_background_data(fits_data, padding, title)
    logging.debug(f"plotting {[len(x) for x in star_lists]} stars per color")
    positions = []
    for stars in star_lists:
        stars = add_pixels(stars, wcs, PADDING)
        positions.append([(o.xpos, o.ypos) for o in stars])
    from itertools import cycle

    cycol = cycle("rgbcmyk")

    for idx, pos in enumerate(positions):
        if len(pos) > 0:
            apps = CircularAperture(pos, r=sizes[idx])
            apps.plot(color=next(cycol), lw=1.5, alpha=0.5)
    random.seed(42)

    if annotate:
        for idx, star in enumerate(star_lists):
            annotate_it(star, -10, 15, random_offset=random_offset[idx], size=10)
    return fig


#  plot_fits is false if no background needs to be plotted, in that case all zeros are used as data
def get_plot_with_background_data(fits_data: List[float], padding: int, title: str):
    fig = plt.figure(figsize=(36, 32), dpi=80, facecolor="w", edgecolor="k")
    plt.title(title, fontsize=40)
    median = np.median(fits_data)
    fits_data = np.pad(
        fits_data, (padding, padding), "constant", constant_values=(100, 100)
    )
    plt.imshow(
        fits_data, cmap="gray_r", origin="lower", vmin=0, vmax=min(median * 5, 65536)
    )
    return fig, fits_data


def save(fig, path):
    fig.savefig(path)
    plt.close(fig)


def run_standard_field_charts(
    star_descriptions: StarDescriptionList,
    wcs,
    fieldchartsdirs,
    reference_fits_frame,
    comp_stars: ComparisonStars,
):
    # setting the font size for titles/axes
    plt.rcParams.update({"axes.titlesize": "large", "axes.labelsize": "large"})
    fits_data, _, _ = reading.get_fits_data(reference_fits_frame)
    fits_data_blank, _, _ = reading.get_fits_data(reference_fits_frame, blank_data=True)
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
    all_stars_no_label = set_custom_label(all_stars_descr, "")

    # field chart with all detections
    logging.info("Plotting field chart with all detected stars...")
    fig = plot_it(
        [all_stars_no_label],
        [4.0],
        [False],
        fits_data,
        wcs,
        "All detected stars",
        PADDING,
        annotate=False,
    )
    save(
        fig, fieldchartsdirs + "all_detections_{}_stars".format(len(all_stars_no_label))
    )

    # vsx stars get their aavso id label
    vsx_descr = utils.get_stars_with_metadata(star_descriptions, "VSX")
    vsx_labeled = set_aavso_id_label(vsx_descr)

    # field chart with all vsx stars
    logging.info("Plotting field chart with all VSX variable stars...")
    fig = plot_it(
        [vsx_labeled], [10.0], [False], fits_data, wcs, "All VSX stars", PADDING
    )
    save(fig, fieldchartsdirs + "vsx_stars_{}".format(len(vsx_labeled)))

    # field chart with all vsx stars without the background
    logging.info(
        "Plotting field chart with all VSX variable stars without reference field..."
    )
    fig = plot_it(
        [vsx_labeled],
        [10.0],
        [False],
        fits_data_blank,
        wcs,
        "VSX without background",
        PADDING,
    )
    save(fig, fieldchartsdirs + "vsx_stars_no_ref_{}".format(len(vsx_labeled)))

    # field chart with only the background
    logging.info("Plotting field chart with only the reference field...")
    fig, _ = get_plot_with_background_data(fits_data, 0, "Reference frame")
    save(fig, fieldchartsdirs + "only_ref")

    # candidate stars get their local id label
    candidate_descr = utils.get_stars_with_metadata(
        star_descriptions, "CANDIDATE", exclude=["VSX"]
    )
    candidate_labeled = set_local_id_label(candidate_descr)

    # localid/radec inputfile stars get their local id label
    selected_desc = utils.get_stars_with_metadata(star_descriptions, "SITE")
    selected_no_vsx_descr = utils.get_stars_with_metadata(
        star_descriptions, "SITE", exclude=["VSX"]
    )
    selected_no_vsx_labeled = set_custom_label(
        selected_no_vsx_descr,
        [x.get_metadata("SITE").our_name for x in selected_no_vsx_descr],
    )
    selected_count = len(selected_no_vsx_labeled)

    # field chart with all vsx stars + candidates + radeccatalog
    logging.info("Plotting field chart with all VSX variable stars + candidate vars...")
    fig = plot_it(
        [vsx_labeled, candidate_labeled],
        [10.0, 5.0],
        [False, True],
        fits_data,
        wcs,
        "VSX stars + candidate stars",
        PADDING,
    )
    save(
        fig,
        fieldchartsdirs
        + f"vsx_{len(vsx_labeled)}_and_candidates_{len(candidate_labeled)}",
    )

    # field chart with all vsx stars + localid/radec inputfile
    logging.info(
        f"Plotting field chart with all VSX variable stars + {selected_count} selected vars..."
    )
    fig = plot_it(
        [vsx_labeled, selected_no_vsx_labeled],
        [10.0, 5.0],
        [False, True],
        fits_data,
        wcs,
        "VSX stars + selected stars",
        PADDING,
    )
    save(
        fig,
        fieldchartsdirs
        + "vsx_{}_and_selected_{}".format(len(vsx_labeled), selected_count),
    )

    # compstars get their vmag as label
    # comp_stars_descr = comp_stars.star_descriptions
    # comp_stars_labeled = set_custom_label(comp_stars_descr, [x.vmag for x in comp_stars_descr])

    logging.info(
        f"Plotting field chart for each of the {selected_count} selected stars"
    )
    # Plotting finder charts for the site
    # field charts for each individually selected starfile star
    for star in tqdm.tqdm(selected_desc, desc="Field chart for each star"):
        filtered_compstars, check_star = do_compstars.filter_comparison_stars(
            star, comp_stars
        )
        filtered_compstars_sds = filtered_compstars.star_descriptions
        check_star_sd = check_star.star_descriptions
        compstars_labeled = set_custom_label(
            filtered_compstars_sds,
            [x.get_metadata("UCAC4").vmag for x in filtered_compstars_sds],
        )
        checkstar_labeled = set_custom_label(
            check_star_sd, f"Kmag={check_star_sd[0].get_metadata('UCAC4').vmag}"
        )
        starui: utils.StarUI = utils.get_star_or_catalog_name(star, suffix="")
        fig = plot_it(
            [[star], vsx_labeled, compstars_labeled, checkstar_labeled],
            [7.0, 5.0, 3.0, 4.0],
            [True, False, True, False],
            fits_data,
            wcs,
            f"VSX stars + comp stars + {starui.catalog_name} (star {star.local_id})",
            PADDING,
        )
        save(fig, f"{fieldchartsdirs}vsx_and_star_{starui.filename_no_ext}")
        gc.collect()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="munipack automation field charts")
    parser.add_argument(
        "-d",
        "--datadir",
        help="The directory where the data can be found (fits in ./fits dir under the data dir",
        nargs="?",
        required=True,
    )
    parser.add_argument("-s", "--stars", help="List the star id's to plot", nargs="+")
    parser.add_argument("-n", "--novsx", help="Don't plot vsx stars", nargs="+")
    args = parser.parse_args()

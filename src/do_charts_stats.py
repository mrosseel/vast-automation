import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import pandas as pd
import glob
from tqdm import tqdm
import logging
import multiprocessing as mp
import reading
from pathlib import Path
from typing import List, Dict
import utils
import do_aavso_report
from star_description import StarDescription
from star_metadata import CompStarData
import re
import toml
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
from utils import StarDict

""" Create charts showing statistics on the detected stars, variables, ... """


def get_fig_and_ax():
    fig = plt.figure(figsize=(20, 12), dpi=150)
    # plt.rcParams['figure.constrained_layout.use'] = True
    plt.gcf().subplots_adjust(bottom=0.2)
    ax = plt.subplot(111)
    return fig, ax


def plot_comparison_stars(
    chartsdir: str,
    stars: List[StarDescription],
    stardict: StarDict,
    jdfilter: List[float],
):
    # for every selected star make one chart
    for star in tqdm(stars, desc="Plotting comparison stars", unit="star"):
        compstars: CompStarData = star.get_metadata("COMPSTARS")
        compstar_ids = compstars.compstar_ids + [compstars.extra_id, star.local_id]
        labels = compstars.compstar_ids + ["K"]
        dfs = reading.read_lightcurve_ids(compstar_ids, stardict)
        # compstars without the 'K' star or check star
        star.result["compA"] = helper_plot_stars(
            star,
            chartsdir,
            utils.get_star_or_catalog_name(star, suffix="_compstarsA"),
            dfs[:-1],
            labels,
            jdfilter,
            show_error=True,
        )
        # compstars + 'K' star
        star.result["compB"] = helper_plot_stars(
            star,
            chartsdir,
            utils.get_star_or_catalog_name(star, suffix="_compstarsB"),
            dfs,
            labels + ["V"],
            jdfilter,
            show_error=False,
        )


def helper_plot_stars(
    star: StarDescription,
    chartsdir: str,
    starui: utils.StarUI,
    dfs: List[pd.DataFrame],
    labels: List[str],
    jdfilter: List[float],
    show_error: bool = False,
):
    """ helper function to plot the comparison stars  """
    fig, ax = get_fig_and_ax()
    xmax, xmin, ymax, ymin = float("-inf"), float("inf"), float("-inf"), float("inf")
    cmap = plt.get_cmap("Set1")
    number = len(dfs)
    colors = [cmap(i) for i in np.linspace(0, 1, number + 1)]
    title = f"{starui.filename_no_suff_no_ext} comp. stars{' + V' if not show_error else ''}"
    for idx, dfx in enumerate(dfs):
        dfx["floatJD"] = dfx["JD"].astype(float)
        dfx["Vrel30"] = dfx["Vrel"] + 30
        df = utils.jd_filter_df(dfx, jdfilter)
        xmax = max(xmax, df["floatJD"].max())
        xmin = min(xmin, df["floatJD"].min())
        ymax = max(ymax, df["Vrel30"].max())
        ymin = min(ymin, df["Vrel30"].min())
        if show_error:
            ax.errorbar(
                df["floatJD"],
                df["Vrel30"],
                yerr=df["err"],
                linestyle="",
                color=colors[idx],
                ms=2,
            )
        else:
            fmt = "*" if labels[idx] == "K" else "^" if labels[idx] == "V" else "."
            ax.plot(df["floatJD"], df["Vrel30"], fmt, color=colors[idx], markersize=2)
    fontp = FontProperties()
    fontp.set_size("18")
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    ax.legend(
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.4),
        fancybox=True,
        shadow=True,
        ncol=len(labels),
        prop=fontp,
    )
    ax.set_title(title)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin - 0.05, ymax + 0.05)
    plt.xlabel("JD")
    plt.ylabel("Instr. mag")
    plt.xticks(rotation=25)
    ax.invert_yaxis()
    plt.autoscale()
    save_location = Path(chartsdir, starui.filename_no_ext + ".png")
    fig.savefig(save_location)
    plt.close(fig)
    return save_location


# read all image.cat.info files
# get JD and aperture via regex
# plot
def plot_aperture_vs_jd(chartsdir: str, vastdir: str, jdfilter: List[float]):
    x, y = get_aperture_and_jd(vastdir, jdfilter)
    logging.info(f"Aperture lengths after filtering: {len(x)} {len(y)}, with jdfilter: {jdfilter}")
    fig, ax = get_fig_and_ax()
    ax.plot(x, y, "*r", markersize=2)
    ax.set_title("Aperture vs JD")
    plt.xlabel("JD")
    plt.ylabel("Aperture (px)")
    plt.xticks(rotation=25)
    plt.minorticks_on()
    plt.autoscale()
    save_location = Path(chartsdir, "aperture_vs_jd" + ".png")
    fig.savefig(save_location)
    plt.close(fig)


def get_catinfo_files(vastdir):
    return [Path(f) for f in glob.glob(vastdir + "/*.cat.info")]


def get_aperture_and_jd(vastdir: str, jdfilter: List[float]):
    catinfo_files = get_catinfo_files(vastdir)
    x = []
    y = []
    for file in tqdm(catinfo_files, desc="Reading apertures per image", unit="image"):
        a, b = get_jd_aperture_from_catinfo(file)
        x.append(float(a))
        y.append(float(b))
    logging.debug(f"Aperture lengths before filtering: {len(x)} {len(y)}")
    return utils.jd_filter_array(x, y, jdfilter)


def plot_merr_vs_jd(chartsdir: str, stars: List[StarDescription], jdfilter):
    dfs = reading.read_lightcurve_sds(stars)
    for star, df in tqdm(
        zip(stars, dfs),
        desc="Plotting magnitude error vs jd",
        unit="star",
        total=len(stars),
    ):
        df["floatJD"] = df["JD"].astype(np.float)
        df = utils.jd_filter_df(df, jdfilter)
        starui: utils.StarUI = utils.get_star_or_catalog_name(star)
        fig, ax = get_fig_and_ax()
        ax.plot(df["floatJD"], df["err"], "*r", markersize=2)
        ax.set_title("Magnitude error vs JD")
        plt.xlabel("JD (day)")
        plt.ylabel("Mag error (mag)")
        plt.xticks(rotation=25)
        plt.minorticks_on()
        plt.autoscale()
        save_location = Path(chartsdir, f"{starui.filename_no_ext}_merr_vs_jd" + ".png")
        fig.savefig(save_location)
        plt.close(fig)
        star.result["merr_vs_jd"] = save_location


# part of plot_apertures
def get_jd_aperture_from_catinfo(filename):
    the_regex = re.compile(
        r"^write_string_to_log_file\(\): JD=\s*(.*)\s*ap=\s*(\S*)\s*.*$"
    )
    with open(filename, "r") as infile:
        for line in infile:
            thesearch = the_regex.search(line)
            if thesearch:
                return thesearch.group(1), thesearch.group(2)
    return None


def plot_aperture_vs_airmass(chartsdir: str, vastdir: str, wcs, jdfilter: List[float]):
    x, y = get_aperture_and_jd(vastdir, jdfilter)
    settings = toml.load("settings.txt")
    sitelat = settings["sitelat"]
    sitelong = settings["sitelong"]
    sitealt = settings["sitealt"]
    earth_location = EarthLocation(lat=sitelat, lon=sitelong, height=sitealt * u.m)

    xair = []
    central_coord = SkyCoord.from_pixel(
        wcs.pixel_shape[0] / 2.0, wcs.pixel_shape[1] / 2.0, wcs=wcs, origin=0
    )
    for entry in tqdm(x, total=len(x), desc="Calculating airmass", unit="JD"):
        result = do_aavso_report.calculate_airmass(
            central_coord, earth_location, entry
        ).value
        xair.append(result)

    fig, ax = get_fig_and_ax()
    ax.plot(xair, y, "*r", markersize=2)
    ax.tick_params(axis="x", which="minor", bottom=False)
    ax.set_title("Aperture vs Airmass")
    plt.xlabel(f"Airmass of central point: {utils.get_hms_dms(central_coord)}")
    plt.ylabel("Aperture (px)")
    plt.minorticks_on()
    plt.autoscale()
    save_location = Path(chartsdir, "aperture_vs_airmass" + ".png")
    fig.savefig(save_location)
    plt.close(fig)


def plot_fwhm(fwhm):
    fig_size = (36, 32)
    dpi = 100
    matplotlib.rcParams.update({"font.size": 38})

    fig, ax = plt.subplots(figsize=fig_size, dpi=dpi, facecolor="w", edgecolor="k")
    ax.set_xlabel("Image number")
    ax.set_ylabel("FWHM")
    ax.set_title(r"Progress over time of the FWHM per image")
    ax.grid(True)
    xaxis = range(0, len(fwhm))
    plt.bar(x=xaxis, height=np.take(fwhm, 1, axis=1), width=1)
    plt.show()
    save(fig, settings.fieldchartsdirs + "fwhm.png")


def save(fig, path):
    fig.savefig(path)
    plt.close(fig)


def main(fwhm):
    plot_fwhm(fwhm)

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


def plot_comparison_stars(chartsdir: str, stars: List[StarDescription], stardict: StarDict):
    # for every selected star make one chart
    for star in tqdm(stars, desc="Plotting comparison stars", unit="star"):
        compstars: CompStarData = star.get_metadata('COMPSTARS')
        compstar_ids = compstars.compstar_ids + [compstars.extra_id, star.local_id]
        labels = compstars.compstar_ids + ['K']
        dfs = reading.read_lightcurve_ids(compstar_ids, stardict)
        star.result['compA'] = \
            plot_star_fluctuations(star, chartsdir, utils.get_star_or_catalog_name(star, suffix="_compstarsA"),
                                   dfs[:-1], labels, show_error=True)
        star.result['compB'] = \
            plot_star_fluctuations(star, chartsdir,  utils.get_star_or_catalog_name(star, suffix="_compstarsB"), dfs,
                                   labels + ['V'], show_error=False)


def plot_star_fluctuations(star: StarDescription, chartsdir: str, starui: utils.StarUI, dfs: List[pd.DataFrame],
                           labels: List[str], show_error: bool = False):
    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = plt.subplot(111)
    xmax, xmin, ymax, ymin = float('-inf'), float('inf'), float('-inf'), float('inf')
    cmap = plt.get_cmap('Set1')
    number = len(dfs)
    colors = [cmap(i) for i in np.linspace(0, 1, number)]
    title = f"{starui.filename_no_suff_no_ext} comp. stars{' + V' if not show_error else ''}"
    for idx, df in enumerate(dfs):
        df['floatJD'] = df['JD'].astype(float).to_numpy()
        df['Vrel30'] = df['Vrel'] + 30
        xmax = max(xmax, df['floatJD'].max())
        xmin = min(xmin, df['floatJD'].min())
        ymax = max(ymax, df['Vrel30'].max())
        ymin = min(ymin, df['Vrel30'].min())
        if show_error:
            ax.errorbar(df['floatJD'], df['Vrel30'], yerr=df['err'], linestyle='', color=colors[idx], ms=2)
        else:
            fmt = '*r' if labels[idx] == 'K' else '^b' if labels[idx] == 'V' else '.'
            ax.plot(df['floatJD'], df['Vrel30'], fmt, color=colors[idx], markersize=2)
    fontp = FontProperties()
    fontp.set_size('18')
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])
    # Put a legend below current axis
    ax.legend(labels, loc='upper center', bbox_to_anchor=(0.5, -0.18), fancybox=True, shadow=True, ncol=len(labels),
              prop=fontp)
    ax.set_title(title)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin - 0.05, ymax + 0.05)
    plt.xlabel('JD')
    plt.ylabel('Instr. mag')
    ax.invert_yaxis()
    save_location = Path(chartsdir, starui.filename_no_ext + ".png")
    fig.savefig(save_location)
    plt.close(fig)
    return save_location


# read all image.cat.info files
# get JD and aperture via regex
# plot
def plot_aperture_vs_jd(chartsdir: str, vastdir: str):
    x, y = get_aperture_and_jd(vastdir)
    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = plt.subplot(111)
    ax.plot(x, y, '*r', markersize=2)
    ax.set_title('Aperture vs JD')
    plt.xlabel('JD')
    plt.ylabel('Aperture (px)')
    plt.minorticks_on()
    save_location = Path(chartsdir, 'aperture_vs_jd' + '.png')
    fig.savefig(save_location)
    plt.close(fig)


def get_aperture_and_jd(vastdir: str):
    catinfo_files = get_catinfo_files(vastdir)
    x = []
    y = []
    for file in tqdm(catinfo_files, desc="Reading apertures per image", unit="image"):
        a, b = get_jd_aperture_from_catinfo(file)
        x.append(float(a))
        y.append(float(b))
    return x, y


def plot_merr_vs_jd(chartsdir: str, stars: List[StarDescription]):
    dfs = reading.read_lightcurve_sds(stars)
    for star, df in tqdm(zip(stars, dfs), desc="Plotting magnitude error vs jd", unit="star", total=len(stars)):
        df['floatJD'] = df['JD'].astype(np.float)
        starui: utils.StarUI = utils.get_star_or_catalog_name(star)
        fig = plt.figure(figsize=(20, 12), dpi=150)
        ax = plt.subplot(111)
        ax.plot(df['floatJD'], df['err'], '*r', markersize=2)
        ax.set_title('Magnitude error vs JD')
        plt.xlabel('JD (day)')
        plt.ylabel('Mag error (mag)')
        plt.minorticks_on()
        save_location = Path(chartsdir, f'{starui.filename_no_ext}_merr_vs_jd' + '.png')
        fig.savefig(save_location)
        plt.close(fig)
        star.result['merr_vs_jd'] = save_location


# part of plot_apertures
def get_catinfo_files(vastdir):
    return [Path(f) for f in glob.glob(vastdir + '/*.cat.info')]


# part of plot_apertures
def get_jd_aperture_from_catinfo(filename):
    the_regex = re.compile(r'^write_string_to_log_file\(\): JD=\s*(.*)\s*ap=\s*(\S*)\s*.*$')
    catalog_dict = {}
    with open(filename, 'r') as infile:
        for line in infile:
            thesearch = the_regex.search(line)
            if thesearch:
                return thesearch.group(1), thesearch.group(2)
    return None


def plot_aperture_vs_airmass(chartsdir: str, vastdir: str, wcs):
    x, y = get_aperture_and_jd(vastdir)
    settings = toml.load('settings.txt')
    sitelat = settings['sitelat']
    sitelong = settings['sitelong']
    sitealt = settings['sitealt']
    earth_location = EarthLocation(lat=sitelat, lon=sitelong, height=sitealt * u.m)

    xair = []
    central_coord = SkyCoord.from_pixel(wcs.pixel_shape[0] / 2.0, wcs.pixel_shape[1] / 2.0, wcs=wcs, origin=0)
    for entry in tqdm(x, total=len(x), desc="Calculating airmass", unit="JD"):
        result = do_aavso_report.calculate_airmass(central_coord, earth_location, entry).value
        xair.append(result)

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = plt.subplot(111)
    ax.plot(xair, y, '*r', markersize=2)
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_title('Aperture vs Airmass')
    plt.xlabel(f'Airmass of central point: {utils.get_hms_dms(central_coord)}')
    plt.ylabel('Aperture (px)')
    plt.minorticks_on()
    save_location = Path(chartsdir, 'aperture_vs_airmass' + '.png')
    fig.savefig(save_location)
    plt.close(fig)


def plot_cumul_histo_detections(savefig=True):
    result = reading.read_lightcurve_ids()
    matplotlib.rcParams.update({'font.size': 38})
    fig_size = (38, 32)
    dpi = 100
    # print(len(result))
    keys = result.keys()
    values = list(map(lambda x: x[0] / x[1] * 100, result.values()))
    # print(len(keys), len(values))
    num_bins = 10
    fig, ax = plt.subplots(figsize=fig_size, dpi=dpi, facecolor='w', edgecolor='k')

    # the histogram of the data
    n, bins, patches = ax.hist(values, num_bins, density=1)

    # add a 'best fit' line
    # y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
    #     np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
    # ax.plot(bins, range(1,50), '--')
    ax.set_xlabel('% of star detections, on which star was found')
    ax.set_ylabel('Nr of stars')
    ax.set_title(r'Cumulative histogram of star detections')
    ax.grid(True)
    # plt.xticks(np.arange(0, 110, step=10))
    major_ticks = np.arange(0, 11000, 2000)
    minor_ticks = np.arange(0, 11000, 1000)
    ax.set_yticks(major_ticks)
    ax.set_yticks(minor_ticks, minor=True)
    plt.minorticks_on()
    plt.xlim(0, 100)
    plt.ylim(0, 10000)
    # plt.yticks(np.arange(0, 11000, step=1000))
    # Tweak spacing to prevent clipping of ylabel
    # fig.tight_layout()
    plt.hist(bins=num_bins, x=values, cumulative=1)
    plt.show()
    if savefig:
        save(fig, settings.fieldchartsdirs + 'cumul_histo_detections.png')

    fig, ax = plt.subplots(figsize=fig_size, dpi=dpi, facecolor='w', edgecolor='k')
    ax.set_xlabel('star number')
    ax.set_ylabel('Star is in % of images')
    ax.set_title(r'Barchart of on which % of images the star is seen')
    ax.grid(True)
    xaxis = range(0, len(values))
    logging.info(f"len xaxis: {len(xaxis)}, len values: {len(values)}")
    plt.bar(x=xaxis, height=sorted(values), width=1)
    plt.show()
    save(fig, settings.fieldchartsdirs + 'barcharts.png')


def plot_fwhm(fwhm):
    fig_size = (36, 32)
    dpi = 100
    matplotlib.rcParams.update({'font.size': 38})

    fig, ax = plt.subplots(figsize=fig_size, dpi=dpi, facecolor='w', edgecolor='k')
    ax.set_xlabel('Image number')
    ax.set_ylabel('FWHM')
    ax.set_title(r'Progress over time of the FWHM per image')
    ax.grid(True)
    xaxis = range(0, len(fwhm))
    plt.bar(x=xaxis, height=np.take(fwhm, 1, axis=1), width=1)
    plt.show()
    save(fig, settings.fieldchartsdirs + 'fwhm.png')


def save(fig, path):
    fig.savefig(path)
    plt.close(fig)


def main(fwhm):
    plot_cumul_histo_detections()
    plot_fwhm(fwhm)

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
from star_description import StarDescription
from star_metadata import CompStarData

StarDict = Dict[int, StarDescription]


#### Create charts showing statistics on the detected stars, variables, ...


def plot_comparison_stars(chartsdir: str, stars: List[StarDescription], stardict):
    # for every selected star make one chart
    for star in stars:
        _, _, _, filename_no_ext = utils.get_star_or_catalog_name(star, suffix='_compstars')
        compstars: CompStarData = star.get_metadata('COMPSTARS')
        compstar_ids = compstars.compstar_ids + [compstars.extra_id]
        labels = compstars.compstar_ids + ['K']
        dfs = read_lightcurves(compstar_ids, stardict)
        plot_star_fluctuations(chartsdir, filename_no_ext, dfs, labels, use_mean=False)
        plot_star_fluctuations(chartsdir, filename_no_ext, dfs, labels, use_mean=True)


def plot_star_fluctuations(chartsdir: str, filename_no_ext: str, dfs, labels: List[str], use_mean: bool = False):
    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = plt.subplot(111)
    xmax, xmin, ymax, ymin = float('-inf'), float('inf'), float('-inf'), float('inf')
    for df in dfs:
        df['JDF'] = df['JD'].astype(float).to_numpy()
        # df['JDF'] = np.arange(0, len(df['JD']))
        df['Vrelrel'] = df['Vrel'] - df['Vrel'].mean() if use_mean else df['Vrel'] + 30
        xmax = max(xmax, df['JDF'].max())
        xmin = min(xmin, df['JDF'].min())
        ymax = max(ymax, df['Vrelrel'].max())
        ymin = min(ymin, df['Vrelrel'].min())
        ax.errorbar(df['JDF'], df['Vrelrel'], yerr=df['err'], linestyle='none')
        # ax.plot(df['JDF'], df['Vrelrel'])
    fontp = FontProperties()
    fontp.set_size('20')
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])
    # Put a legend below current axis
    print(f"labels size is {len(labels)}")
    ax.legend(labels, loc='upper center', bbox_to_anchor=(0.5, -0.18), fancybox=True, shadow=True, ncol=len(labels),
              prop=fontp)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin - 0.5, ymax + 0.5)
    plt.xlabel('JD')
    plt.ylabel('Instr. mag')
    ax.invert_yaxis()
    save_location = Path(chartsdir, filename_no_ext + f"{'A' if use_mean else 'B'}.png")
    print(f"writing {save_location}")
    fig.savefig(save_location)
    plt.close(fig)


def plot_cumul_histo_detections(savefig=True):
    result = read_lightcurves()
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


def old_read_lightcurves(lightcurvedir: str, nr_threads):
    files = glob.glob(lightcurvedir + '*.txt')
    result = {}
    pool = mp.Pool(nr_threads * 2, maxtasksperchild=None)

    for file, partial_result in tqdm(pool.imap_unordered(read_lightcurve, files), total=len(files),
                                     desc='Reading all lightcurves'):
        result[file] = partial_result

    return result


def read_lightcurves(star_ids: List[int], stardict: StarDict):
    result = []
    for star_id in star_ids:
        df = reading.read_lightcurve_vast(stardict[star_id].path)
        result.append(df)
    return result


def read_lightcurve(file):
    df = pd.read_csv(file, skiprows=[1], sep=' ')
    length = len(df.index)
    df = df[df['V-C'] < 99]
    filterlength = len(df.index)
    # print (filterlength/length)
    return file, [filterlength, length]


def save(fig, path):
    fig.savefig(path)
    plt.close(fig)


def main(fwhm):
    plot_cumul_histo_detections()
    plot_fwhm(fwhm)

import pandas as pd
import numpy as np
import init_loader
from init_loader import init, settings, meta_init
import reading
import logging
import glob
import argparse
import os
# from importlib import reload
# main_muniwin = reload(main_muniwin)
import do_calibration
from reading import file_selector
from reading import trash_and_recreate_dir
from reading import create_dir
from do_aperture import gather_data

# globals photometry
maxstar = -1
apertures = "2,2.73,3.82,5.27,7.09,9.27,11.82,14.73,18,21.64,25.64,30"
gain = 1.4  # ADC gain of Josch's camera
skyinner = 8
skyouter = 14

# globals matching
max_stars = 30
vertices = 6
clip_thresh = 1000
sp_fields = 0
sp_maxoffset = 2000


# find_best_x => perform_x ->

# get a subset of the files
# generate X number of photometry config files
# generate phot files for each of the config files
# get non-infinite counts
# print results
def find_best_photometry_config(percentage):
    # get subset of files for photometry
    trash_and_recreate_dir(settings.testdir)
    pattern = 'kout*.fts'
    logging.info(f"Selecting files from {settings.convfitsdir} with pattern {pattern}")
    chosen_files = file_selector(settings.convfitsdir, pattern, percentage)
    logging.info(f"Chosen {len(chosen_files)} files per photometry run, files: {chosen_files}")
    filelist = generate_photometry_config_files()
    perform_photometry(filelist, chosen_files)


#
def find_best_matching_config(winning_photometry_dir, percentage, resume=False):
    all_phot_files = glob.glob(winning_photometry_dir + "*.pht")
    logging.info(f"Chosen {len(all_phot_files)} photometry files from dir {winning_photometry_dir}")
    filelist = generate_matching_config_files()
    reference_frame = all_phot_files[0]
    perform_matching(filelist, all_phot_files, reference_frame)


# generate the config files
def write_photometry_config(filename, apertures, maxstar, fwhm, thresh, gain, skyinner, skyouter):
    a = {'apertures': apertures, 'maxstar': maxstar, 'fwhm': fwhm, 'thresh': thresh, 'gain': gain, 'skyinner': skyinner,
         'skyouter': skyouter}
    write_config(filename, a)


# generate the config files
def write_matching_config(filename, max_stars, vertices, clip_thresh, sp_fields, sp_maxoffset):
    a = {'max_stars': max_stars, 'vertices': vertices, 'clip_thresh': clip_thresh, 'sp_fields': sp_fields,
         'sp_maxoffset': sp_maxoffset}
    write_config(filename, a)


def write_config(filename, dict_of_entries):
    l = []
    for k, v in dict_of_entries.items():
        l.append((k, v))
    b = pd.DataFrame(l)
    b.to_csv(filename, header=False, index=False, sep="=", mode='a')


def generate_photometry_config_files():
    logging.info("Generating photometry config files...")
    create_dir(settings.testdir + 'conf')
    result = []
    counter = 0
    for fwhm in np.arange(1, 8, 0.5):
        for thresh in range(1, 8, 1):
            counter += 1
            filename = f'{settings.testdir}conf/muniphot{counter}.conf'
            write_photometry_config(filename, apertures, maxstar, fwhm, thresh, gain, skyinner, skyouter)
            result.append(filename)
    logging.info(f"{len(result)} config files generated.")
    return result


# max_stars=30
# vertices=6
# clip_thresh=1000
# sp_fields=0
# sp_maxoffset=2000
def generate_matching_config_files():
    logging.info("Generating matching config files...")
    trash_and_recreate_dir(settings.testdir + 'matching')
    result = []
    counter = 0
    for max_stars in np.arange(10, 60, 15):
        for vertices in range(3, 10, 2):
            for clip_thresh in range(100, 5000, 1000):
                for sp_fields in range(0, 3, 1):
                    for sp_maxoffset in range(1, 3000, 1000):
                        counter += 1
                        filename = f'{settings.testdir}matching/match{counter}.conf'
                        logging.info(f"Writing to filename {filename}")
                        write_matching_config(filename, max_stars, vertices, clip_thresh, sp_fields, sp_maxoffset)
                        result.append(filename)
    logging.info(f"{len(result)} matching config files generated.")
    return result


# offset param is to allow generating non-standard photometry in a range beyond the existing dirs
def perform_photometry(configfilelist, chosen_files, offset=1):
    # filestring = ' '.join(chosen_files)
    logging.debug(f"Perform photometry {chosen_files}, {configfilelist}")
    for id, entry in enumerate(configfilelist):
        logging.debug(f"Perform photometry {id} {entry}")
        outputdir = f'{settings.testdir}{id + offset:05d}/'
        create_dir(outputdir)
        main_muniwin.write_photometry(config_file=entry, files=chosen_files, outputdir=outputdir)


# offset param is to allow generating non-standard photometry in a range beyond the existing dirs
def perform_matching(configfilelist, all_photometry_files, reference_frame, offset=1):
    input_file_string = ' '.join(all_photometry_files)
    # either read the previous reference frame or calculate a new one
    logging.debug(f"Perform photometry {input_file_string}, {configfilelist}")
    for id, entry in enumerate(configfilelist):
        logging.debug(f"Perform photometry {id} {entry}")
        outputdir = f'{settings.testdir}match{id + offset:05d}/'
        create_dir(outputdir)
        logging.info("========================================")
        for photfile in all_photometry_files:
            main_muniwin.write_match(photfile, reference_frame, config_file=entry, to_match_is_full_path=True,
                                     outputdir=outputdir)


def analyse(resultdirs=None, apertureidx=None):
    logging.info("starting analysis")
    if resultdirs is None:
        logging.info("resultdirs is none")
        resultdirs = sorted([f for f in glob.glob(settings.testdir + "*" + os.path.sep) if not "conf" in f])
        logging.info(resultdirs)
    if apertureidx is None:
        _, apertureidx, _ = reading.read_aperture()
    logging.info(f"apertureidx: {apertureidx}")
    df = pd.DataFrame(columns=['resultdir', 'nrstars', 'starsdetectedpct', 'realpercentage'])
    for resultdir in resultdirs:
        # reading photometry files
        logging.info(f"Reading dir: {resultdir}")
        jdphot, fwhmphot, nrstarsphot, star_result = read_photometry.read_photometry(star_list_1=init.star_list,
                                                                                     apertureidx=apertureidx,
                                                                                     matched_files=glob.glob(
                                                                                         resultdir + '/*.pht'),
                                                                                     fake_references=True)
        logging.info(
            f"Results: nr stars photometry:{len(nrstarsphot)}, sum of {np.sum(nrstarsphot)}, length of star list:{len(init.star_list)}")
        starsdetectedpct = np.sum(nrstarsphot) / (len(init.star_list) * len(nrstarsphot)) * 100
        logging.info(star_result.flatten())
        logging.info(f"number of finite star results: {np.sum(np.isfinite(star_result.flatten().shape))}")
        logging.info(f"Percent of stars detected: {starsdetectedpct}")
        df.loc[len(df)] = [resultdir, np.sum(nrstarsphot.flatten()), starsdetectedpct,
                           show_percentage_real(star_result)]
    df.to_csv(settings.testdir + 'results_phot.txt', index=False)
    winningdf = df.sort_values(['realpercentage', 'nrstars'], ascending=[False, False]).iloc[0]
    logging.info(f"Best result: {winningdf}")
    return winningdf['resultdir']


def analyse_match(resultdirs=None, apertureidx=None):
    if resultdirs is None:
        resultdirs = sorted([f for f in glob.glob(settings.testdir + "match*" + os.path.sep) if not "matching" in f])
        logging.info(resultdirs)
    if apertureidx is None:
        _, apertureidx, _ = reading.read_aperture()
    logging.info(f"apertureidx: {apertureidx}")
    df = pd.DataFrame(columns=['resultdir', 'nrstars', 'starsdetectedpct', 'realpercentage'])
    for resultdir in resultdirs:
        # reading photometry files
        logging.info(f"Reading dir:  {resultdir}")
        jdphot, fwhmphot, nrstarsphot, star_result = read_photometry.read_photometry(star_list_1=init.star_list,
                                                                                     apertureidx=apertureidx,
                                                                                     matched_files=glob.glob(
                                                                                         resultdir + '/*.pht'),
                                                                                     fake_references=False)
        logging.info(f"Results: nr stars photometry:{len(nrstarsphot)}, " +
                     f"sum of {np.sum(nrstarsphot)}, length of star list:{len(init.star_list)}")
        starsdetectedpct = np.sum(nrstarsphot) / (len(init.star_list) * len(nrstarsphot)) * 100
        logging.info(star_result.flatten())
        logging.info(f"number of finite star results: {np.sum(np.isfinite(star_result.flatten().shape))}")
        logging.info(f"Percent of stars detected: {starsdetectedpct}")
        df.loc[len(df)] = [resultdir, np.sum(nrstarsphot.flatten()), starsdetectedpct,
                           show_percentage_real(star_result)]
    df.to_csv(settings.testdir + 'results_match.txt', index=False)
    winningdf = df.sort_values(['realpercentage', 'nrstars'], ascending=[False, False]).iloc[0]
    logging.info(f"Best result: {winningdf}")
    return winningdf['resultdir']


def show_percentage_real(star_result):
    totalcount = 0
    for fileentry in star_result:
        totalcount += np.sum(np.isfinite(fileentry)) / 2
    result = totalcount / (len(init.star_list) * len(star_result))
    logging.info(f"Percentage real: {result}")
    return result


# ax = s.hist()  # s is an instance of Series
# fig = ax.get_figure()
# fig.savefig('/path/to/figure.pdf')

# takes a star_list and a dir, and returns a reduced star list - all stars which already have a file in that dir are removed
# def reduce_matching_config_list(config_list, the_path):
#     the_dir = os.listdir(the_path)
#     the_dir.sort()
#     found = []
#     for filename in the_dir:
#         found.append(filename_to_star(filename))
#     print("Found", len(found), "stars already processed in", the_path)
#     return [item for item in star_list_1 if item not in found]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Testing photometry and matching config params')
    parser.add_argument('percentage')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (fits in ./fits dir under the data dir",
                        nargs='?', required=True)
    args = parser.parse_args()
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    logging.info("before init")
    meta_init(args.datadir)
    # global init
    init = init_loader.init
    settings = init_loader.settings
    import main_muniwin
    import read_photometry

    # from importlib import reload
    # main_muniwin = reload(main_muniwin)

    maxstar = len(init.star_list)
    find_best_photometry_config(args.percentage)
    resultdir = analyse(apertureidx=None)
    # find_best_matching_config(resultdir, args.percentage, resume=True)
    # analyse_match(apertureidx=init.aperture)

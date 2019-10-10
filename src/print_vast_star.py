import argparse
import logging
import main_vast
import do_calibration
import utils

def main(args):
    vastdir = utils.add_trailing_slash(args.datadir)
    wcs_file = vastdir+'new-image.fits'
    wcs = do_calibration.get_wcs(wcs_file)
    all_files = main_vast.file_selector(the_dir=vastdir, match_pattern="*.dat")
    file_targets = [vastdir+star_to_dat(int(x)) for x in args.stars]
    print("file targets ", file_targets)
    # selected_files = [x for x in all_files where]
    all_stardict = main_vast.read_stardict(vastdir)
    logging.info(f"Number of found lightcurves: {len(file_targets)}, number of identified stars: {len(all_stardict.keys())}")
    args.upsilon = args.vsx = args.selectedstarfile = False
    sds = main_vast.construct_star_descriptions(vastdir, wcs, all_stardict, file_targets, args)
    for star in sds:
        print(star, '\n')

def star_to_dat(star: int):
    return f"out{star:05}.dat"

if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description='Print all data on one star')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (usually the vast dir)",
                        nargs='?', required=True)
    parser.add_argument('-s', '--stars', help="List the star id's to plot", nargs='+')
    parser.add_argument('-p', '--plot', help="Plot the stars using pgfv")
    args = parser.parse_args()
    main(args)

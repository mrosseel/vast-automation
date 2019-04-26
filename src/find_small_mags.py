import do_aperture as dop
from init_loader import init, settings
import glob
import argparse
import logging
import time
from tqdm import tqdm

def process_one_pht(the_file, apertureidx:int, star0=None):
    logging.debug("in process")
    start = time.perf_counter()
    photheader, apertures, nrstars, stars, stardata = dop.read_pht_file('', open(the_file, 'rb').read(),
                                                                        only_apertureidx=apertureidx, read_stars=False)
    end = time.perf_counter()
    result = stardata

    if star0:
        logging.info(f"stardata for star {star0}: {stardata[int(star0)]}")
    else:
        for idx, entry in enumerate(stardata):
            # if entry.mag > 99:
            #     print(f"found big thing: {the_file}")
            if 0.000001 > entry.mag and -0.000001 < entry.mag:
                logging.info(f"found small thing: {the_file} at star0 index:{idx}")
    return end-start

def do_all(apertureidx, directory, star0=None):
    phtdir = directory + "*.pht"
    logging.debug(f"init dir is {settings.matchedphotometrydir}")
    logging.info(phtdir)

    files = sorted(glob.glob(phtdir))
    #print(files)
    total_count = 0
    for entry in tqdm(files, total=len(files)):
        total_count += process_one_pht(entry, apertureidx, star0)
    logging.info(total_count)

def do_one(apertureidx, directory, the_file, star0):
    phtfile= directory + the_file
    logging.debug(f"Loading file at {phtfile}")
    process_one_pht(phtfile, apertureidx, star0)

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    parser = argparse.ArgumentParser(description='Finding strange magnitudes in .pht files')
    parser.add_argument('apertureidx')
    parser.add_argument('directory')
    parser.add_argument('--file')
    parser.add_argument('--star0')
    args = parser.parse_args()
    if args.star0 and not args.file:
        do_all(int(args.apertureidx), args.directory, args.star0)
    if args.file and args.star0:
        do_one(int(args.apertureidx), args.directory, args.file, int(args.star0))
    else:
        do_all(int(args.apertureidx), args.directory)

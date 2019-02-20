import do_photometry as dop
import init
import glob
import argparse
import logging
import time
from tqdm import tqdm

def process_one_pht(the_file, apertureidx: int):
    logging.debug("in process")
    start = time.perf_counter()
    photheader, apertures, nrstars, stars, stardata = dop.read_pht_file('', open(the_file, 'rb').read(),
                                                                        only_apertureidx=apertureidx, read_stars=True)
    end = time.perf_counter()
    result = stardata
    # result = [x[apertureidx] for x in stardata]
    #print(stardata)
    if any(None in x for x in result):
        print("None found in ", the_file)
    if any(0.000001 > x.mag and -0.000001 < x.mag for x in result):
        print(f"found small thing: {the_file}")
    for entry in stardata:
        if entry.mag > 99:
            print(f"found big thing: {the_file}")
        if 0.000001 > entry.mag and -0.000001 < entry.mag:
            print(f"found small thing: {the_file}")

    # a = list(map(lambda x: list(map(lambda y: y.mag, x)), stardata))
    # print(a)
    # b = reduce(lambda x,y: x+y,a)
    # list(filter(lambda x: x<1 or x>99, b))
    return end-start

def main(apertureidx, directory):
    phtdir = directory + "*.pht"
    logging.debug(f"init dir is {init.matchedphotometrydir}")
    print(phtdir)

    files = sorted(glob.glob(phtdir))
    #print(files)
    total_count = 0
    for entry in tqdm(files, total=len(files)):
        total_count += process_one_pht(entry, apertureidx)
    print(total_count)

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    parser = argparse.ArgumentParser(description='Finding strange magnitudes in .pht files')
    parser.add_argument('apertureidx')
    parser.add_argument('directory')
    args = parser.parse_args()
    main(int(args.apertureidx), args.directory)

import do_photometry as dop
import init
import glob
import argparse
import logging
import timeit
from tqdm import tqdm

def process_one_pht(the_file, apertureidx: int):
    logging.debug("in process")
    # start = timeit.timeit()
    photheader, apertures, nrstars, stars, stardata = dop.read_pht_file('', open(the_file, 'rb').read(),
                                                                        only_apertureidx=apertureidx)
    # end = timeit.timeit()
    result = [x[apertureidx] for x in stardata]
    # print(end - start)
    # print(result)
    if any(None in x for x in result):
        print("None found in ", the_file)
    if any(0.000001 > x.mag and -0.000001 < x.mag for x in result):
        print(f"found small thing: {the_file}")
    # a = list(map(lambda x: list(map(lambda y: y.mag, x)), stardata))
    # print(a)
    # b = reduce(lambda x,y: x+y,a)
    # list(filter(lambda x: x<1 or x>99, b))


def main(apertureidx):
    phtdir = init.matchedphotometrydir + "*.pht"
    logging.debug(f"init dir is {init.matchedphotometrydir}")
    print(phtdir)

    files = sorted(glob.glob(phtdir))
    #print(files)
    for entry in tqdm(files, total=len(files)):
        process_one_pht(entry, apertureidx)

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    parser = argparse.ArgumentParser(description='Finding strange magnitudes in .pht files')
    parser.add_argument('apertureidx')
    args = parser.parse_args()
    main(int(args.apertureidx))


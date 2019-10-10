from astropy.io import fits
from functools import partial
import glob
import logging
import multiprocessing as mp
import tqdm

# TODO all C's need to be changed into CV

# read the filter value for each fits file, to be used in aavso reporting

def read_filters(files=None, inputdir: str = None):
    if inputdir is None:
        inputdir = settings.convfitsdir
    if files is None:
        files = glob.glob(inputdir + "*.fts")
    logging.info(f"Reading filters for {len(files)} files...")
    pool = mp.Pool(init.nr_threads * 2, maxtasksperchild=100)
    func = partial(read_fits_header)
    result = {}
    logging.debug(f"files: {files}")
    for filter_tuple in tqdm.tqdm(pool.imap_unordered(func, files, 5), total=len(files)):
        result[filter_tuple[0]] = filter_tuple[1]
    logging.debug(f"read_filters result: {result}")
    return result

def read_fits_header(file):
    fts = fits.open(file)
    jd = fts[0].header['JD']
    result_filter = fts[0].header['FILTER']
    return jd, result_filter.strip()


if __name__ == '__main__':
    print("not implemented yet")
    # init_loader.meta_init('./current/')
    # # global init
    # init = init_loader.init
    # settings = init_loader.settings
    # logger = logging.getLogger()
    # logger.setLevel(logging.DEBUG)
    # logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    # result = read_filters()


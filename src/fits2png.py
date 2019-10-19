import astropy_helper
import logging
import argparse
import utils
from pathlib import PurePath, Path
from datetime import datetime
import glob
import subprocess


def process_fits(data_dir: Path, result_dir: Path, start: str, end: str, extension='*.fit'):
    selected_fits = data_dir.glob(extension)
    # print("first 10 fits:", list(selected_fits)[:10])
    result = []
    for fits in selected_fits:
        _, _, _, _, jd = astropy_helper.get_data_from_fits_header(open(fits, 'rb'))
        result.append((float(jd), fits))
    result = sorted(result, key=lambda x: x[0])
    within_date_range = list(filter(lambda x: start <= x[0] <= end, result))
    for entry in within_date_range:
        print("processing", entry)
        subprocess.run(["fits2bitmap", "--stretch", "log", entry[1], "-o",
                        resultdir / f"{str(entry[0])}_{entry[1].with_suffix('.png').name}"])


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description='munipack automation cli')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (usually the vast dir)",
                        nargs='?', required=True)
    parser.add_argument('-r', '--resultdir',
                        help="The directory where all results will be written",
                        nargs='?', required=True)
    parser.add_argument('-s', '--start',
                        help="The start Julian Date",
                        nargs='?', required=True)
    parser.add_argument('-e', '--end',
                        help="The end Julian Date",
                        nargs='?', required=True)
    parser.add_argument('-x', '--verbose', help="Set logging to debug mode", action="store_true")
    args = parser.parse_args()
    datadir = Path(args.datadir)
    datenow = datetime.now()
    resultdir = Path(args.resultdir)
    filehandler = f"{datadir}vastlog-{datenow:%Y%M%d-%H_%M_%S}.log"
    fh = logging.FileHandler(filehandler)
    fh.setLevel(logging.INFO)
    # add the handlers to the logger
    logger.addHandler(fh)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        fh.setLevel(logging.DEBUG)

    process_fits(datadir, resultdir, float(args.start), float(args.end))

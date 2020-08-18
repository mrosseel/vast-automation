from multiprocessing import cpu_count
import logging
import argparse
import utils
from datetime import datetime
import main_vast
import os.path

if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description='munipack automation cli')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (usually the vast dir)",
                        required=True)
    parser.add_argument('-r', '--resultdir',
                        help="The directory where all results will be written",
                        required=True)
    parser.add_argument('-k', '--checkstarfile', help="The bright and stable stars used to do ensemble photometry",
                        required=False)
    parser.add_argument('-o', '--owncatalog', help="Supply a file to identify stars known to you by RA/DEC",
                        required=False)
    parser.add_argument('-s', '--selectedstarfile',
                        help="Load a file with star ids, these ids will be used for field charts/reporting")
    parser.add_argument('--apikey', help="Astrometry.net api key for automatic plate solving of the referenc "
                                               "frame found by vast.",
                        required=False)
    parser.add_argument('--fitsdir', help="The dir where the fits are, only needed to plate-solve the reference frame",
                        required=False)
    parser.add_argument('-v', '--vsx', help="Add vsx stars to field charts/reporting list", action="store_true")
    parser.add_argument('-!', '--allstars', help="Generate phase/light/aavso for all stars. WARNING: takes long",
                        action="store_true")
    parser.add_argument('-f', '--field', help="Generate field charts (finder charts)", action="store_true")
    parser.add_argument('--stats', help="Generate stats charts", action="store_true")
    parser.add_argument('-p', '--phase', help="Generate phase charts", action="store_true")
    parser.add_argument('-i', '--light', help="Generate lightcurve charts", action="store_true")
    parser.add_argument('-c', '--candidates', help="Generate phase diagrams for autocandidates", action="store_true")
    parser.add_argument('-a', '--aavso', help="Generate aavso reports",
                        action="store_true")
    parser.add_argument('-t', '--aavsolimit', help="Limits the number of lines per aavso file. -t 5000 splits the"
                                                   "observations in files of 5000 lines each.",
                        nargs='?', type=int, default=None, const=5000)
    parser.add_argument('--threads', help="Specifies the number of threads to be used by this program", nargs='?', type=int, default=cpu_count()-1)
    parser.add_argument('--selectvsx', help="Add all VSX stars to the selected stars", action="store_true")
    parser.add_argument('--selectcandidates', help="Add all candidate stars to the selected stars", action="store_true")
    parser.add_argument('--site', help="Generate a hugo compatible page")
    parser.add_argument('-x', '--verbose', help="Set logging to debug mode", action="store_true")
    parser.add_argument('-l', '--laststars', help="Use the star descriptions of the previous run to do the charting",
                        action="store_true")
    parser.add_argument('-u', '--upsilon', help="Add upsilon star info to charting", action="store_true")
    parser.add_argument('--jdfilter', help="Filters out a range of JD's. Filters every JD between the first and "
                        "the second argument. If you want to filter from a certain point up to infinity,"
                        "use a suitably large JD like 9999999", nargs='+', type=float, required=False)
    parser.add_argument('--jdrefignore', help="Igone that the ref frame jd is inside of the filter", action="store_true")
    args = parser.parse_args()
    datadir = utils.add_trailing_slash(args.datadir)
    datenow = datetime.now()
    filehandler = f"{datadir}vastlog-{datenow:%Y%M%d-%H_%M_%S}.log"
    fh = logging.FileHandler(filehandler)
    fh.setLevel(logging.INFO)
    # add the handlers to the logger
    logger.addHandler(fh)

    # sanity checks
    # improve with: https://stackoverflow.com/questions/37471636/python-argument-parsing-validation-best-practices
    assert os.path.exists(args.datadir), "datadir does not exist"
    # assert os.path.exists(args.resultdir), "resultdir does not exist" ==> this dir is created
    assert os.path.exists(args.checkstarfile) if args.checkstarfile else True, "checkstarfile does not exist"
    assert os.path.exists(args.owncatalog) if args.owncatalog else True, "owncatalog does not exist"
    assert os.path.exists(args.selectedstarfile) if args.selectedstarfile else True, "selectedstarfile does not exist"
    if args.selectcandidates and not args.candidates:
        logging.warning("--selectcandidates was selected, but --candidates was absent. Corrected.")
        args.candidates = True

    if args.selectvsx and not args.vsx:
        logging.warning("--selectvsx was selected, but --vsx was absent. Corrected.")
        args.vsx = True

    # monitoring_thread = start_monitoring(seconds_frozen=15, test_interval=1000)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        fh.setLevel(logging.DEBUG)

    logging.info(f"Do phase: {'YES' if args.phase else 'NO'}, Do light: {'YES' if args.light else 'NO'},"
                 f"Do field: {'YES' if args.field else 'NO'}, Do aavso: {'YES' if args.aavso else 'NO'},"
                 f"Do site: {'YES' if args.site else 'NO'}")
    main_vast.run_do_rest(args)

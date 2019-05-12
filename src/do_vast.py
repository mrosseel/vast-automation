import init_loader
from init_loader import init, settings
import logging
import argparse
import subprocess
import utils
from datetime import datetime

if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description='munipack automation cli')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (fits in ./fits dir under the data dir",
                        nargs='?', required=True)
    parser.add_argument('-n', '--nowait', help="Don't wait 10 secs before starting", action="store_true")
    parser.add_argument('-v', '--vsx', help="Add vsx stars to field charts/reporting list", action="store_true")
    parser.add_argument('-s', '--starfile', help="Load a file with star ids, these ids will be used for field charts/reporting")
    parser.add_argument('-x', '--verbose', help="Set logging to debug mode", action="store_true")
    parser.add_argument('-l', '--laststars', help="Use the star descriptions of the previous run to do the charting",
                        action="store_true")
    parser.add_argument('-u', '--upsilon', help="Add upsilon star info to charting", action="store_true")
    args = parser.parse_args()
    datadir = utils.add_trailing_slash(args.datadir)
    datenow = datetime.now()
    filehandler = f"{datadir}vastlog-{datenow:%Y%M%d-%H_%M_%S}.log"
    fh = logging.FileHandler(filehandler)
    fh.setLevel(logging.INFO)
    # add the handlers to the logger
    logger.addHandler(fh)
    init_loader.meta_init(datadir)
    # global init
    init = init_loader.init
    settings = init_loader.settings
    # print(dir(init))
    # print(dir(settings))pr
    import main_vast

    #monitoring_thread = start_monitoring(seconds_frozen=15, test_interval=1000)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        fh.setLevel(logging.DEBUG)

    logging.info(f"Calculating {len(init.star_list)} stars from base dir: {settings.basedir} \
          \nupsilon:\t{init.do_ml} \
          \nlight plot:\t{init.do_lightcurve_plot} \
          \nphasediagram:\t{init.do_phase_diagram} \
          \nfield charts:\t{init.do_field_charts} \
          \nreporting:\t{init.do_reporting}")
    if not args.nowait:
        logging.info("Press Enter to continue...")
        subprocess.call("read -t 10", shell=True, executable='/bin/bash')
    main_vast.run_do_rest(init.do_ml, init.do_lightcurve_plot, init.do_phase_diagram,
                             init.do_field_charts, init.do_reporting, args)

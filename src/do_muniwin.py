import init_loader
from init_loader import init, settings
import logging
import do_lightcurve
import argparse
import subprocess
import utils
from datetime import datetime
#from hanging_threads import start_monitoring

if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description='munipack automation cli')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (fits in ./fits dir under the data dir",
                        nargs='?', required=True)
    parser.add_argument('-c', '--chart', help="Generate lightcurve, lightcurve plot and phase diagram plot", nargs='+')
    parser.add_argument('-n', '--nowait', help="Don't wait 10 secs before starting", action="store_true")
    parser.add_argument('-v', '--vsx', help="Only do charting for the vsx stars", action="store_true")
    parser.add_argument('-l', '--laststars', help="Use the star descriptions of the previous run to do the charting",
                        action="store_true")
    parser.add_argument('-u', '--upsilon', help="Add upsilon star info to the star descriptions", action="store_true")
    args = parser.parse_args()
    datadir = utils.add_trailing_slash(args.datadir)

    fh = logging.FileHandler(f"{datadir}munilog-{datetime.now():%Y%M%d-%H_%M_%S}.log")
    fh.setLevel(logging.INFO)
    # add the handlers to the logger
    logger.addHandler(fh)
    init_loader.meta_init(datadir)
    # global init
    init = init_loader.init
    settings = init_loader.settings
    # print(dir(init))
    # print(dir(settings))pr
    import main_muniwin

    #monitoring_thread = start_monitoring(seconds_frozen=15, test_interval=1000)

    if args.chart:
        logging.info("in chart part")
        logging.info(f"Writing lightcurves... {[x.local_id for x in star_descriptions_ucac4]}")
        chosen_stars = [x.local_id for x in star_descriptions_ucac4]
        do_lightcurve.write_lightcurves(chosen_stars,
                                  comparison_stars_1, aperture, int(apertureidx), jd, fwhm, star_result)
        logging.info("starting charting / phase diagrams...")
        logging.info(f"comparison stars decs: {comparison_stars_1_desc}")
        do_charts.run(star_descriptions_ucac4, comparison_stars_1_desc, do_lightcurve_plot, do_phase_diagram)

    else:
        logging.info(f"Calculating {len(init.star_list)} stars from base dir: {settings.basedir} \
              \nconvert_fits:\t{init.do_convert_fits} \
              \nphotometry:\t{init.do_photometry} \
              \nmatch:\t\t{init.do_match} \
              \naperture:\t{init.do_aperture_search} \
              \npos:\t\t{init.do_pos} \
              \ncopstars:\t{init.do_compstars} \
              \nlightcurve:\t{init.do_lightcurve} \
              \nupsilon:\t{init.do_ml} \
              \nlight plot:\t{init.do_lightcurve_plot} \
              \nphasediagram:\t{init.do_phase_diagram} \
              \nfield charts:\t{init.do_field_charts} \
              \nreporting:\t{init.do_reporting}")
        if not args.nowait:
            logging.info("Press Enter to continue...")
            subprocess.call("read -t 10", shell=True, executable='/bin/bash')
        main_muniwin.run_do_rest(init.do_convert_fits, init.do_photometry, init.do_match, init.do_compstars,
                                 init.do_aperture_search, init.do_lightcurve,
                                 init.do_pos, init.do_ml, init.do_lightcurve_plot, init.do_phase_diagram,
                                 init.do_field_charts, init.do_reporting, args)

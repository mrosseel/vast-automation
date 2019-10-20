import aavso
import reading
import logging
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
# supress the warning about vector transforms so as not to clutter the doc build log
import warnings

warnings.filterwarnings('ignore', module='astropy.coordinates.baseframe')
from star_description import StarDescription
import argparse
import read_camera_filters
from typing import Tuple, List
from pandas import DataFrame
from comparison_stars import ComparisonStars
from pathlib import PurePath


def calculate_airmass(coord, location, jd):
    time = Time(jd, format='jd')
    alt_azs = coord.transform_to(AltAz(obstime=time, location=location))
    return alt_azs.secz


def report(star: StarDescription, df_curve: DataFrame, target_dir: PurePath, vastdir: str, sitelat, sitelong,
           sitealt, comparison_stars: ComparisonStars, filter=None, observer='RMH', chunk_size=None):
    comparison_star = comparison_stars.star_descriptions[0]
    star_match_ucac4, separation = star.get_match_string("UCAC4")
    star_match_vsx, separation = star.get_match_string("VSX", strict=False)
    comp_ucac4 = comparison_star.get_match_string("UCAC4", strict=True)
    var_display_name = star_match_ucac4 if star_match_vsx is None else star_match_vsx
    var_display_name = var_display_name if var_display_name is not None else f"Star_{star.local_id}"
    check_display_name = comparison_star.aavso_id if comparison_star.aavso_id is not None else comp_ucac4[0]

    # logging.info(" Star match:{}, comparison_star:{}".format(var_display_name, comparison_star))
    comparison_star_vmag = comparison_star.vmag
    # NO COMP STAR TODAY
    comparison_star_vmag = 0.0
    title = f"{star.local_id:05}" if star.aavso_id is None else star.aavso_id
    earth_location = EarthLocation(lat=sitelat, lon=sitelong, height=sitealt * u.m)
    logging.debug(f"Starting aavso report with star:{star}")
    if chunk_size is None:
        chunk_size = df_curve.shape[0]
    star_chunks = [df_curve[i:i + chunk_size] for i in range(0, df_curve.shape[0], chunk_size)]
    chunk_counters = 0

    filterdict = None
    # Setting up the filter value
    if filter is None:
        filterdict = read_camera_filters.read_filters()
        print(filterdict)
        filterlambda = lambda x: filterdict[x]
    else:
        filterlambda = lambda x: filter

    for chunk in star_chunks:
        chunk_counters += 1
        suffix = f"_{chunk_counters}.txt" if len(star_chunks) != 1 else ".txt"
        with open(PurePath(target_dir, f"{title}_ext{suffix}"), 'w') as fp:
            writer = aavso.ExtendedFormatWriter(fp, observer, software='munipack-automation', type='EXTENDED',
                                                obstype='CCD')
            for _, row in chunk.iterrows():
                # logging.info(row, type(row))
                writer.writerow({
                    'name': var_display_name,
                    'date': row['JD'],
                    'magnitude': row['realV'],
                    'magnitude_error': row['realErr'],
                    'filter': filterlambda(row['JD']),
                    'transformed': 'YES',
                    'magnitude_type': 'STD',
                    'comparison_name': 'ENSEMBLE',
                    'comparison_magnitude': 'na',
                    'check_name': check_display_name,
                    'check_magnitude': comparison_star_vmag,
                    'airmass': calculate_airmass(star.coords, earth_location, row['floatJD']),
                    'group': 'na',
                    'chart': 'na',
                    'notes': 'na'
                })


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description='munipack automation cli')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (fits in ./fits dir under the data dir",
                        nargs='?', required=True)
    parser.add_argument('-c', '--chart', help="Only generate lightcurve, lightcurve plot and phase diagram plot",
                        nargs='+')
    parser.add_argument('-n', '--nowait', help="Don't wait 10 secs before starting", action="store_true")
    parser.add_argument('-v', '--vsx', help="Add vsx stars to field charts/reporting list", action="store_true")
    parser.add_argument('-s', '--starfile',
                        help="Load a file with star ids, these ids will be used for field charts/reporting")
    parser.add_argument('-x', '--verbose', help="Set logging to debug mode", action="store_true")
    parser.add_argument('-l', '--laststars', help="Use the star descriptions of the previous run to do the charting",
                        action="store_true")
    parser.add_argument('-u', '--upsilon', help="Add upsilon star info to charting", action="store_true")
    args = parser.parse_args()
    datadir = utils.add_trailing_slash(args.datadir)
    datenow = datetime.now()
    filehandler = f"{datadir}munilog-{datenow:%Y%M%d-%H_%M_%S}.log"
    fh = logging.FileHandler(filehandler)
    fh.setLevel(logging.INFO)
    # add the handlers to the logger
    logger.addHandler(fh)
    # print(dir(init))
    # print(dir(settings))pr
    import main_muniwin

    # monitoring_thread = start_monitoring(seconds_frozen=15, test_interval=1000)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        fh.setLevel(logging.DEBUG)

    # star_descriptions_ucac4 = do_calibration.add_ucac4_to_star_descriptions(star_descriptions)
    # logging.info(f"AAVSO Reporting with: {len(selected_stars)} stars")
    # trash_and_recreate_dir(settings.aavsoreportsdir)
    # for star in selected_stars:
    #     do_aavso_report.report(settings.aavsoreportsdir, star, comparison_stars_1_desc[0])

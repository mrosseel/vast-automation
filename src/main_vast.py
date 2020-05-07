from multiprocessing import cpu_count
import logging
import re
import os
import os.path
import sys
import numpy as np
import time
from collections import namedtuple
import do_calibration
import do_charts_vast
import do_charts_field
import do_charts_stats
import do_compstars
import reading
import utils
import utils_sd
from utils import get_localid_to_sd_dict
from star_description import StarDescription
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.wcs import WCS
from typing import List, Dict, Tuple
from comparison_stars import ComparisonStars
from pathlib import PurePath, Path
from ucac4 import UCAC4
import hugo_site
import pandas as pd
import toml
import subprocess
from star_metadata import CatalogData, SiteData, CompStarData, SelectedFileData
from utils import StarDict

vsx_catalog_name = "vsx_catalog.bin"
vsxcatalogdir = PurePath(os.getcwd(), vsx_catalog_name)
STAR_KEEPER_PERCENTAGE = 0.1


def run_do_rest(args):
    thread_count = cpu_count() - 1
    vastdir = utils.add_trailing_slash(args.datadir)
    resultdir = clean_and_create_resultdir(args.resultdir, vastdir)
    fieldchartsdir = resultdir + 'fieldcharts/'
    do_light = args.light
    do_phase = args.phase
    do_aavso = args.aavso
    logging.info(f"Dir with VaST files: '{vastdir}', results dir: '{resultdir}'")
    # get wcs model from the reference header. Used in writing world positions and field charts
    wcs_file, wcs = reading.read_wcs_file(vastdir)
    ref_jd, _, _, reference_frame = reading.extract_reference_frame(vastdir)
    _, _, _, first_frame = reading.extract_first_frame(vastdir)
    referene_frame_path = Path(reference_frame)
    reference_frame_filename = referene_frame_path.name
    logging.info(f"The reference frame is '{reference_frame}' at JD: {ref_jd}")
    logging.info(f"The first frame is '{first_frame}'")
    logging.info(f"Reference header is '{wcs_file}'")
    # Log filtering settings + check that reference frame is not inside of filter
    if args.jdfilter:
        logging.info(f"Filtering JD's: {args.jdfilter}")
        ref_inside_filter = float(ref_jd) > args.jdfilter[0] and float(ref_jd) > args.jdfilter[1]
        if ref_inside_filter:
            if not args.jdrefignore:
                assert float(ref_jd) > args.jdfilter[0] and float(ref_jd) > args.jdfilter[1], \
                    "Reference frame JD is filtered"
            else:
                logging.info("Reference frame JD is inside of the JD filter")

    if not os.path.isfile(wcs_file):
        full_ref_path = Path(args.fitsdir) / reference_frame_filename
        if not args.fitsdir and args.apikey:
            logging.error("There is no plate-solved reference frame {wcs_file}, please specify both --apikey "
                          "and --fitsdir.")
            sys.exit(0)
        rotation = reading.extract_reference_frame_rotation(vastdir, reference_frame_filename)
        assert rotation == 0.0, f"Error: rotation is {rotation} and should always be 0.0"
        subprocess.Popen(f"python3 ./src/astrometry_api.py --apikey={args.apikey} "
                         f"--upload={full_ref_path} --newfits={wcs_file} --private --no_commercial", shell=True)
        while not os.path.isfile(wcs_file):
            logging.info(f"Waiting for the astrometry.net plate solve...")
            time.sleep(10)

    star_descriptions = construct_star_descriptions(vastdir, resultdir, wcs, args)
    stardict = get_localid_to_sd_dict(star_descriptions)
    logging.debug(f"First (max) 10 star descriptions: "
                  f"{star_descriptions[:10] if (len(star_descriptions) >= 10) else star_descriptions}")
    write_augmented_autocandidates(vastdir, resultdir, stardict)
    write_augmented_all_stars(vastdir, resultdir, stardict)
    owncatalog = utils.get_stars_with_metadata(star_descriptions, "OWNCATALOG")
    logging.info(f"There are {len(owncatalog)} own catalog stars")
    candidate_stars = utils.get_stars_with_metadata(star_descriptions, "CANDIDATE", exclude=["VSX"])
    candidate_stars = utils.add_star_lists(candidate_stars, owncatalog)
    if args.selectcandidates:
        tag_candidates_as_selected(candidate_stars)
    logging.info(f"There are {len(candidate_stars)} candidate stars")

    vsx_stars = utils.get_stars_with_metadata(star_descriptions, "VSX")
    logging.info(f"There are {len(vsx_stars)} vsx stars")
    selected_stars = utils.get_stars_with_metadata(star_descriptions, "SELECTEDFILE")
    # if args.selectvsx:
    #     selected_stars = utils.concat_sd_lists(selected_stars, vsx_stars)
    logging.info(f"There are {len(selected_stars)} selected stars")
    compstar_needing_stars = utils.concat_sd_lists(selected_stars, vsx_stars, candidate_stars, owncatalog)
    comp_stars = set_comp_stars_and_ucac4(star_descriptions, selected_stars, args.checkstarfile, vastdir, stardict,
                                          ref_jd)
    # Set comp stars for all interesting stars (stars which are interesting enough to measure)
    logging.info("Setting per star comparison stars...")
    if args.checkstarfile:
        utils.add_metadata(star_descriptions, CompStarData(compstar_ids=comp_stars.ids))
    else:
        do_compstars.add_closest_compstars(compstar_needing_stars, comp_stars, 10)

    logging.info(f"Using {thread_count} threads for phase plots, lightcurves, ...")
    if args.allstars:
        do_charts_vast.run(star_descriptions, comp_stars, vastdir, resultdir, 'phase_all/', 'light_all/', 'aavso_all/',
                           do_phase=do_phase, do_light=do_light, do_aavso=do_aavso, nr_threads=thread_count,
                           jdfilter=args.jdfilter, desc="Phase/light/aavso of ALL stars")
    else:
        if args.vsx:
            logging.info(f"Plotting {len(vsx_stars)} vsx stars...")
            do_charts_vast.run(vsx_stars, comp_stars, vastdir, resultdir, 'phase_vsx/', 'light_vsx/', 'aavso_vsx/',
                               do_phase=do_phase, do_light=do_light, do_light_raw=do_light,
                               do_aavso=do_aavso, nr_threads=thread_count,
                               jdfilter=args.jdfilter, desc="Phase/light/aavso of VSX stars")
        if args.selectedstarfile or args.owncatalog:
            do_charts_vast.run(selected_stars, comp_stars, vastdir, resultdir, 'phase_selected/', 'light_selected/',
                               'aavso_selected', do_phase=do_phase, do_light=do_light, do_light_raw=do_light,
                               do_aavso=do_aavso, nr_threads=thread_count,
                               jdfilter=args.jdfilter,
                               desc="Phase/light/aavso of selected stars")
        if args.candidates:
            logging.info(f"Plotting {len(candidate_stars)} candidates...")
            do_charts_vast.run(candidate_stars, comp_stars, vastdir, resultdir, 'phase_candidates/',
                               'light_candidates/', 'aavso_candidates/', do_phase=do_phase, do_light=do_light,
                               do_light_raw=do_light, do_aavso=do_aavso, nr_threads=thread_count,
                               jdfilter=args.jdfilter, desc="Phase/light/aavso of candidates")

    # starfiledata is filled in during the phase plotting, so should come after it. Without phase it will be incomplete
    ids = [x.local_id for x in selected_stars]
    logging.info(f"Writing selected files with {len(selected_stars)}  selected stars: {ids}")
    write_selected_files(resultdir, vastdir, selected_stars)
    if args.field:
        do_charts_field.run_standard_field_charts(star_descriptions, wcs, fieldchartsdir, wcs_file, comp_stars)

    if args.stats:
        do_charts_stats.plot_comparison_stars(fieldchartsdir, selected_stars, stardict)
        do_charts_stats.plot_aperture_vs_jd(fieldchartsdir, vastdir)
        do_charts_stats.plot_aperture_vs_airmass(fieldchartsdir, vastdir, wcs)
        do_charts_stats.plot_merr_vs_jd(fieldchartsdir, selected_stars)

    if args.site:
        ids = [x.local_id for x in selected_stars]
        logging.info(f"Creating HTML site with {len(selected_stars)} selected stars: {ids}")
        hugo_site.run(args.site, selected_stars, len(vsx_stars), referene_frame_path, resultdir)


# Either read UCAC4 check stars from a file, or calculate our own comparison stars
def set_comp_stars_and_ucac4(star_descriptions: List[StarDescription], selectedstars: List[StarDescription],
                             checkstarfile: str, vastdir: str, stardict: StarDict, ref_jd: str) -> ComparisonStars:
    ucac4 = UCAC4()
    if checkstarfile:
        # load comparison stars
        checkstars = read_checkstars(checkstarfile)
        ucac4.add_ucac4_to_sd(selectedstars)
        comparison_stars_ids, comparison_stars_1_sds = do_compstars.get_fixed_compstars(star_descriptions, checkstars)
    else:
        ucac4.add_ucac4_to_sd(star_descriptions)
        comparison_stars_ids, comparison_stars_1_sds = do_compstars.get_calculated_compstars(vastdir, stardict,
                                                                                             ref_jd)
    # get all observations for the comparison stars
    comp_observations = []
    for star in comparison_stars_ids:
        comp_magdict = reading.read_magdict_for_star(vastdir, star)
        # logging.info(f"Read comp magdict for {star}: {read_comp_magdict}")
        comp_observations.append(comp_magdict)
    comp_catalogmags = []
    comp_catalogerr = []
    for star in comparison_stars_1_sds:
        star_ucac4_catalog = star.get_metadata("UCAC4")
        comp_catalogmags.append(star_ucac4_catalog.vmag)
        comp_catalogerr.append(star_ucac4_catalog.vmag_err)
    comp_stars = ComparisonStars(comparison_stars_ids, comparison_stars_1_sds, comp_observations, comp_catalogmags,
                                 comp_catalogerr)
    logging.info(
        f"Using {len(comparison_stars_ids)} comparison stars with on average "
        f"{np.array(list(map(len, comp_observations))).mean()} observations")
    return comp_stars


def clean_and_create_resultdir(argsdir: str, vastdir: str):
    resultdir = utils.add_trailing_slash(argsdir) if argsdir is not None else vastdir
    # if resultdir does not exist, create it
    if not os.path.isdir(resultdir):
        logging.info(f"The resultdir '{resultdir}' does not exist, creating it...")
        reading.create_dir(resultdir)
    return resultdir


# quickly test a few xy2sky conversions using our wcs and astropy
def wcs_test_pattern(wcs):
    logging.info("Outputting wcs testing pattern")
    test_pattern = [(0, 0), (500, 500), (136, 985), (-50, 500)]
    for tuple in test_pattern:
        logging.info(f"test pattern: {tuple}")
        result = wcs.all_pix2world(tuple[0], tuple[1], 0, ra_dec_order=True)
        logging.info(f"result: {result[0]}, {result[1]}")


def read_checkstars(checkstar_file: str) -> List[str]:
    result = []
    with open(checkstar_file) as file:
        for line in file:
            result.append(line.strip())
    return result


def get_autocandidates(vastdir: str) -> List[int]:
    origname = 'vast_autocandidates.log'
    result = []
    with open(Path(vastdir, origname), 'r', encoding='utf-8') as infile:
        for line in infile:
            linetext = line.rstrip()
            star_id = utils.get_starid_from_outfile(linetext)
            result.append(star_id)
    return result


def write_augmented_autocandidates(readdir: str, writedir: str, stardict: StarDict):
    origname = f"{readdir}vast_autocandidates.log"
    newname = f"{writedir}vast_autocandidates_pos.txt"
    logging.info(f"Writing {newname}...")
    with open(origname, 'r', encoding='utf-8') as infile, open(newname, 'w') as outfile:
        for line in infile:
            linetext = line.rstrip()
            star_id = utils.get_starid_from_outfile(linetext)
            if star_id in stardict:
                cacheentry = stardict[star_id]
                outfile.write(
                    f"{linetext}{'' if cacheentry.path is not '' else '*'}\t{cacheentry.aavso_id}\t{utils.get_lesve_coords(cacheentry.coords)}\n")
            else:
                outfile.write(f"{linetext}*\t{'None'}\n")


def write_augmented_all_stars(readdir: str, writedir: str, stardict: StarDict):
    origname = f"{readdir}vast_list_of_all_stars.log"
    newname = f"{writedir}vast_list_of_all_stars_pos.txt"
    logging.info(f"Writing {newname}...")
    with open(origname, 'r', encoding='utf-8') as infile, open(newname, 'w') as outfile:
        for line in infile:
            star_id = line.split()[0]
            if int(star_id) in stardict:
                cacheentry = stardict[int(star_id)]
                outfile.write(
                    f"{star_id}\t{cacheentry.aavso_id}\t{utils.get_hms_dms(cacheentry.coords)}\t{cacheentry.coords.ra} {cacheentry.coords.dec}\n")


# naam, ra, dec, max, min, type, periode, epoch?
def write_selected_files(resultdir: str, vastdir: str, selected_stars: List[StarDescription]):
    owncatalog = f"{resultdir}selected_radec.txt"
    selectedstars = f"{resultdir}selected_localid.txt"
    logging.info(f"Writing {owncatalog} and {selectedstars} with {len(selected_stars)} stars...")
    sorted_stars = utils.sort_selected(selected_stars)
    vsx_stars_len = len(utils.get_stars_with_metadata(selected_stars, "VSX"))
    no_vsx_len = len(utils.get_stars_with_metadata(selected_stars, "SELECTEDFILE", exclude=["VSX"]))
    with open(owncatalog, 'w') as outowncatalog, open(selectedstars, 'w') as outselected:
        preamble = f"# resultdir: {resultdir}, vastdir: {vastdir}, vsx stars: {vsx_stars_len}, " \
                   f"other stars: {no_vsx_len}\n"
        outowncatalog.write(f"{preamble}# our_name,ra,dec,ucac4_name,ucac4_ra,ucac4_dec,minmax,min,max,var_type,"
                            f"period,period_err,epoch\n")
        outselected.write(f"{preamble}# our_name,local_id,minmax,min,max,var_type,period,period_err,epoch\n")


        def format_float_arg(atoml, arg: str, precision):
            if arg is None or arg not in atoml:
                return ''
            if not isinstance(atoml[arg], float):
                return atoml[arg]
            return f"{atoml[arg]:.{precision}f}"


        def format_float_5(atoml, arg: str):
            return format_float_arg(atoml, arg, 5)


        def format_float_1(atoml, arg: str):
            return format_float_arg(atoml, arg, 1)


        def format_string(arg: str, atoml):
            if arg in atoml:
                return atoml[arg]
            return ''


        for star in sorted_stars:
            metadata: SiteData = star.get_metadata("SITE")
            ucac4: CatalogData = utils.get_ucac4_of_sd(star)
            ucac4_name = "" if not ucac4 else f"{ucac4.name}"
            ucac4_coords = "," if not ucac4 else f"{ucac4.coords.ra.deg:.7f},{ucac4.coords.dec.deg:.7f}"
            starui: utils.StarUI = utils.get_star_or_catalog_name(star)
            txt_path = Path(Path(star.result['phase']).parent, "txt", starui.filename_no_ext + '.txt')
            try:
                parsed_toml = toml.load(txt_path)
                outowncatalog.write(
                    f"{metadata.our_name},{star.coords.ra.deg:.7f},{star.coords.dec.deg:.7f},{ucac4_name},{ucac4_coords},"
                    f"{format_string('minmax', parsed_toml)},{format_float_1(parsed_toml, 'min')},"
                    f"{format_float_1(parsed_toml, 'max')},{metadata.var_type},"
                    f"{format_float_5(parsed_toml, 'period')},{format_float_5(parsed_toml, 'period_err')},"
                    f"{format_string('epoch', parsed_toml)}\n")
                outselected.write(
                    f"{metadata.our_name},{star.local_id},"
                    f"{format_string('minmax', parsed_toml)},{format_float_1(parsed_toml, 'min')},"
                    f"{format_float_1(parsed_toml, 'max')},{metadata.var_type},"
                    f"{format_float_5(parsed_toml, 'period')},{format_float_5(parsed_toml, 'period_err')},"
                    f"{format_string('epoch', parsed_toml)}\n")
            except FileNotFoundError:
                logging.error(f"While writing selected files, Could not find {txt_path}")


# write file with vsx local id, name and ra/dec, and write file to be used as 'selected star' file
def write_vsx_stars(resultdir, results_ids, stars: List[StarDescription]):
    newname = f"{resultdir}vsx_stars.txt"
    selected_file = f"{resultdir}vsx_stars_selected.txt"
    logging.info(f"Writing {newname}...")
    total_found = 0
    stardict = utils.get_localid_to_sd_dict(stars)
    logging.debug(f"Receiving {len(stardict.keys())} as vsx input")
    with open(newname, 'wt') as fp, open(selected_file, 'wt') as selected:
        for number, vsx_id in enumerate(results_ids):
            current_sd = stardict[vsx_id]
            found = False if current_sd.path is '' else True
            assert vsx_id == current_sd.local_id
            total_found += 1 if found else 0
            fp.write(
                f"{vsx_id}{'' if found else '*'}:\t{current_sd.aavso_id}\t{utils.get_lesve_coords(current_sd.coords)}\n")
            selected.write(f'{vsx_id},,VSX-{number},,\n')
        fp.write(
            f"# Total entries: {len(results_ids)}, found: {total_found}, not found: {len(results_ids) - total_found}\n")


def write_candidate_stars(resultdir, stars: List[StarDescription]):
    candidates = utils.get_stars_with_metadata(stars, "CANDIDATES", exclude=["VSX"])
    newname = f"{resultdir}candidate_stars.txt"
    logging.info(f"Writing {newname}...")
    total_found = 0
    stardict = utils.get_localid_to_sd_dict(stars)
    logging.debug(f"Receiving {len(stardict.keys())} as vsx input")
    with open(newname, 'wt') as fp:
        for index, current_sd in enumerate(candidates):
            found = False if current_sd.path is '' else True
            total_found += 1 if found else 0
            fp.write(f'{current_sd.local_id},,CANDIDATE-{index},,\n')


def count_dat_entries(afile):
    return sum(1 for line in open(afile, 'r') if line.rstrip())


def construct_star_descriptions(vastdir: str, resultdir: str, wcs: WCS, args):
    star_descriptions = utils_sd.construct_raw_star_descriptions(vastdir, wcs, None, STAR_KEEPER_PERCENTAGE)
    logging.info(f"Number of stars on more than {STAR_KEEPER_PERCENTAGE:.0%} of frames: {len(star_descriptions)}")
    stardict = get_localid_to_sd_dict(star_descriptions)

    # Add VSX information to SDs
    star_descriptions, results_ids = do_calibration.add_vsx_names_to_star_descriptions(star_descriptions,
                                                                                       vsxcatalogdir, 0.01)
    logging.debug(f"Identified {len(results_ids)} VSX stars")
    vsx_stars = utils.get_stars_with_metadata(star_descriptions, "VSX")
    assert len(vsx_stars) == len(results_ids)
    logging.debug(f"Test Tagged {len(vsx_stars)} stars as VSX.")

    # write the vsx stars used into a file
    results_ids.sort()
    write_vsx_stars(resultdir, results_ids, star_descriptions)

    # tag all candidates with a 'candidate' catalog
    tag_candidates(vastdir, star_descriptions)
    write_candidate_stars(resultdir, star_descriptions)

    # adds sitedata to vsx stars
    if args.selectvsx:
        tag_vsx_as_selected(vsx_stars)

    # adds sitedata to selected stars
    if args.selectedstarfile:
        tag_selected(args.selectedstarfile, stardict)
        logging.debug(
            f"Succesfully read {len(list(filter(lambda x: x.has_metadata('SELECTEDFILE'), star_descriptions)))} "
            f"stars from file:"
            f" {[x.local_id for x in list(filter(lambda x: x.has_metadata('SELECTEDFILE'), star_descriptions))]}")

    if args.owncatalog:
        tag_owncatalog(args.owncatalog, star_descriptions)

    return star_descriptions


def tag_candidates(vastdir: str, star_descriptions: List[StarDescription]):
    candidate_ids = get_autocandidates(vastdir)
    candidate_stars = do_calibration.select_star_descriptions(candidate_ids, star_descriptions)
    do_calibration.add_metadata_to_star_descriptions(candidate_stars, ["CANDIDATE"], strict=False)


def tag_selected(selectedstarfile: str, stardict: StarDict):
    try:
        df = pd.read_csv(selectedstarfile, delimiter=',', comment='#',
                         names=['our_name', 'local_id', 'minmax', 'min', 'max', 'var_type', 'period', 'period_err',
                                'epoch'],
                         dtype={'local_id': int, 'minmax': str, 'epoch': float},
                         skipinitialspace=True, warn_bad_lines=True)
        df = df.replace({np.nan: None})
        logging.info(f"Selecting {len(df)} stars added by {selectedstarfile}: {df['local_id'].to_numpy()}")
        for idx, row in df.iterrows():
            the_star: StarDescription = stardict.get(row['local_id'])
            if the_star is None:
                logging.error(f"Could not find star {row['local_id']}, consider removing it from your txt file")
                continue
            the_star.metadata = SiteData(minmax=row['minmax'],
                                         var_min=row['min'],
                                         var_max=row['max'],
                                         var_type=row['var_type'],
                                         our_name=row['our_name'],
                                         period=float(row['period'])
                                         if row['period'] is not None and row['period'] is not 'None' else None,
                                         period_err=float(row['period_err'])
                                         if row['period_err'] is not None and row['period_err'] is not 'None' else None,
                                         source="OWN",
                                         epoch=row['epoch'])
            the_star.metadata = SelectedFileData()
            logging.debug(f"starfile {the_star.local_id} metadata: {the_star.metadata}, "
                          f"{the_star.get_metadata('SELECTEDFILE')}")
            logging.debug(f"starfile {the_star.get_metadata('SELECTEDFILE')}")
        logging.debug(f"Tagged {len(df)} stars as selected by file.")
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        import traceback
        print(traceback.print_exc())
        logging.error(message)
        logging.error(f"Could not read {selectedstarfile}, star {row['local_id']}")


def tag_vsx_as_selected(vsx_stars: List[StarDescription]):
    for the_star in vsx_stars:
        if the_star.has_metadata("SITE"):  # don't overwrite the SITE entry of SELECTEDFILE which has priority
            continue
        vsx_metadata = the_star.get_metadata("VSX")
        extradata = vsx_metadata.extradata
        if extradata is None:
            logging.error(f"Could not find extradata for star {the_star.local_id}, "
                          f"consider removing it from your txt file")
            continue
        # extradata: {'id': index, 'OID': row['OID'], 'Name': row['Name'], 'Type': row['Type'],
        # 'l_Period': row['l_Period'], 'Period': row['Period'], 'u_Period': row['u_Period']})
        the_star.metadata = SiteData(var_type=str(extradata['Type']),
                                     vsx_var_flag=str(extradata['V']),
                                     separation=float(vsx_metadata.separation),
                                     our_name=str(extradata['Name']),
                                     period=float(extradata['Period'])
                                     if not np.isnan(extradata['Period']) else None,
                                     period_err=extradata['u_Period']
                                     if not extradata['u_Period'] is None or not np.isnan(extradata['u_Period'])
                                     else None,
                                     var_min=float(extradata['min']) if not np.isnan(extradata['min']) else None,
                                     var_max=float(extradata['max']) if not np.isnan(extradata['max']) else None,
                                     minmax=construct_vsx_mag_range(extradata),
                                     source='VSX'
                                     )
        the_star.metadata = SelectedFileData()
        logging.debug(f"site {the_star.local_id} metadata: {the_star.metadata}, "
                      f"{the_star.get_metadata('SITE')}")
        logging.debug(f"site {the_star.get_metadata('SITE')}")
    logging.debug(f"Tagged {len(vsx_stars)} stars as selected vxs stars.")


def tag_candidates_as_selected(candidate_stars: List[StarDescription]):
    for the_star in candidate_stars:
        if the_star.has_metadata("SITE"):  # don't overwrite the SITE entry of SELECTEDFILE or VSX which has priority
            continue
        candidate_metadata = the_star.get_metadata("CANDIDATE")
        if candidate_metadata is None:
            logging.error(f"Could not find CANDIDATE metadata for star {the_star.local_id}")
            continue
        the_star.metadata = SiteData(our_name=str(the_star.local_id), source='CANDIDATE')
        the_star.metadata = SelectedFileData()
        logging.debug(f"site {the_star.local_id} metadata: {the_star.metadata}, "
                      f"{the_star.get_metadata('SITE')}")
        logging.debug(f"site {the_star.get_metadata('SITE')}")
    logging.debug(f"Tagged {len(candidate_stars)} stars as selected candidate stars.")


def construct_vsx_mag_range(entry):
    """ Constructs the min-max rang of a VSX star using an entry from the VSX file """
    def empty_if_nan(x):
        if type(x) is str:
            return x
        return "" if np.isnan(x) else x


    return f"{empty_if_nan(entry['f_min'])} {empty_if_nan(entry['l_min'])} {empty_if_nan(entry['l_max'])} " \
           f"{empty_if_nan(entry['max'])} {empty_if_nan(entry['u_max'])} {empty_if_nan(entry['n_max'])} " \
           f"{empty_if_nan(entry['min'])} {empty_if_nan(entry['u_min'])} {empty_if_nan(entry['n_min'])}"


def tag_owncatalog(owncatalog: str, stars: List[StarDescription]):
    # outfile.write(f"# our_name,ra,dec,minmax,var_type,period,epoch\n")
    logging.info(f"Using owncatalog: {owncatalog}")
    df = pd.read_csv(owncatalog, delimiter=',', comment='#',
                     names=['our_name', 'ra', 'dec', 'ucac4_name', 'ucac4_ra', 'ucac4_dec', 'minmax', 'min', 'max',
                            'var_type', 'period', 'period_err',
                            'epoch'],
                     dtype={'ra': float, 'dec': float, 'ucac4_ra': float, 'ucac4_dec': float, 'minmax': str,
                            'epoch': float},
                     skipinitialspace=True, warn_bad_lines=True)
    df.loc[df['ucac4_ra'].notnull(), 'chosenRA'] = df['ucac4_ra']
    df.loc[df['ucac4_ra'].isnull(), 'chosenRA'] = df['ra']
    df.loc[df['ucac4_dec'].notnull(), 'chosenDEC'] = df['ucac4_dec']
    df.loc[df['ucac4_dec'].isnull(), 'chosenDEC'] = df['dec']
    logging.info(f"tag_owncatalog: {df}")
    ra, dec = (df['chosenRA'], df['chosenDEC'])
    df = df.replace({np.nan: None})
    skycoord: SkyCoord = do_calibration.create_generic_astropy_catalog(ra, dec)
    star_catalog = do_calibration.create_star_descriptions_catalog(stars)
    idx, d2d, d3d = match_coordinates_sky(skycoord, star_catalog, nthneighbor=1)
    for count, index in enumerate(idx):
        row = df.iloc[count]
        the_star = stars[index]
        the_star.metadata = CatalogData(key="OWNCATALOG", catalog_id=row['our_name'],
                                        name=row['our_name'],
                                        coords=SkyCoord(row['ra'], row['dec'], unit="deg"),
                                        separation=d2d[count].degree)
        logging.info(f"Matching {row['our_name']} to {the_star.local_id} with sep {d2d[count].degree} "
                     f"at coords {the_star.coords}")
        if not the_star.has_metadata("SITE"):
            try:
                period = float(row['period']) if row['period'] is not None else None
            except ValueError:
                period = None
            try:
                period_err = float(row['period_err']) if row['period_err'] is not None else None
            except ValueError:
                period_err = None
            try:
                min = float(row['min']) if row['min'] is not None else None
            except ValueError:
                min = None
            try:
                max = float(row['max']) if row['max'] is not None else None
            except ValueError:
                max = None
            the_star.metadata = SiteData(minmax=row['minmax'],
                                         var_min=min,
                                         var_max=max,
                                         var_type=row['var_type'],
                                         our_name=row['our_name'],
                                         period=period,
                                         period_err=period_err,
                                         separation=d2d[count].degree,
                                         source="OWN",
                                         epoch=row['epoch'])
            the_star.metadata = SelectedFileData()
        if d2d[count].degree > 0.01:
            logging.warning(f"Separation between {df.iloc[count]['our_name']} "
                            f"and {stars[index].local_id} is {d2d[count]}")


def set_lines(star: StarDescription):
    with open(star.path) as file:
        star.obs = sum(1 for line in file if line.rstrip())
    return star


def has_option(obj, attr_name):
    return hasattr(obj, attr_name) and obj[attr_name]


def interact():
    import code
    code.InteractiveConsole(locals=dict(globals(), **locals())).interact()

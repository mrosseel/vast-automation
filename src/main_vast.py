import do_calibration
from importlib import reload
import do_charts_vast
import do_charts_field
import reading
reading = reload(reading)
from reading import trash_and_recreate_dir
from star_description import StarDescription
from photometry_blob import PhotometryBlob
import re
import logging
from init_loader import init, settings
from reading import file_selector



def run_do_rest(do_ml, do_lightcurve_plot, do_phase_diagram, do_field_charting, do_reporting, args):

    # either read the previous reference frame or calculate a new one
    # _, _, reference_frame_index = do_calibration.get_reference_frame(100, do_calibration.select_reference_frame_jpeg)

    vast_dir = '/home/jovyan/work/support/vast/vast-1.0rc84/'
    selected_files = file_selector(the_dir=vast_dir, match_pattern="*.dat")
    nr_selected = len(selected_files)

    logging.info(f"Selected {nr_selected} from vast dir.")

    logging.info(f"reference header is {settings.reference_header}")
    # get wcs model from the reference header. Used in writing world positions and field charts
    # wcs = do_calibration.get_wcs(settings.reference_header)
    apertures = None
    aperture = None
    apertureidx = None
    comparison_stars_1, comparison_stars_1_desc = None, None
    photometry_blob = PhotometryBlob() # not yet used


    # construction of the star descriptions list
    star_descriptions = construct_star_descriptions(vast_dir, selected_files, args)

    # select a subset of star_descriptions, as specified in starfile
    # if args.starfile:
    #     selected_filter = partial(catalog_filter, catalog_name='SELECTED')
    #     selected_stars = list(filter(selected_filter, star_descriptions))
    #     logging.debug(f"Light curve: selecting chosen stars with {len(selected_stars)} stars selected.")
    # else:
    #     logging.debug(f"No starfile, selecting all stars")
    #     selected_stars = star_descriptions
    #

    if do_lightcurve_plot or do_phase_diagram:
        logging.info("starting charting / phase diagrams...")
        do_charts_vast.run(star_descriptions, do_lightcurve_plot, do_phase_diagram)

    # if do_field_charting:
    #     logging.info("Starting field chart plotting...")
    #     do_charts_field.run_standard_field_charts(selected_stars, wcs)

    # import code
    # code.InteractiveConsole(locals=dict(globals(), **locals())).interact()
    # if do_reporting:
    #     # star_descriptions_ucac4 = do_calibration.add_ucac4_to_star_descriptions(star_descriptions)
    #     logging.info(f"AAVSO Reporting with: {len(selected_stars)} stars")
    #     trash_and_recreate_dir(settings.aavsoreportsdir)
    #     for star in selected_stars:
    #         do_aavso_report.report(settings.aavsoreportsdir, star, comparison_stars_1_desc[0])


def construct_star_descriptions(main_path, selected_files, args):
    # Start with the list of all detected stars
    star_id_list = []
    selected_files.sort()
    for afile in selected_files:
        m = re.search('out(.*).dat', afile)
        star_id = m.group(1).lstrip('0')
        star_id_list.append(star_id)
    #logging.debug(f"Star id list is : {star_id_list}")

    star_descriptions = do_calibration.get_empty_star_descriptions(star_id_list)
    for idx, sd in enumerate(star_descriptions):
        sd.path = selected_files[idx]

    if args.upsilon:
        # now we generate a list of StarDescriptions, with which all further processing will be done
        logging.info("Setting star_descriptions to upsilon candidates")
        star_descriptions = do_calibration.add_candidates_to_star_descriptions(star_descriptions, 0.1)

    if args.vsx:
        star_descriptions, results_ids = do_calibration.add_vsx_names_to_star_descriptions(star_descriptions, 0.01)
        logging.info(f"Added {len(results_ids)} vsx names to star descriptions")
        do_calibration.add_selected_match_to_stars(star_descriptions, results_ids) # select star ids

    if args.starfile:
        with open(settings.basedir + args.starfile, 'r') as fp:
            lines = fp.readlines()
            starlist = [x.rstrip() for x in lines]
            logging.debug(f"The list of stars read from the starfile is: {starlist} ")
            starlist = [int(x) for x in filter(str.isdigit, starlist)]
            logging.debug(f"The list of stars read from the starfile is: {starlist} ")
            logging.info(f"Selecting {starlist} stars added by {args.starfile}")
            do_calibration.add_selected_match_to_stars(star_descriptions, starlist) # select star ids

    return star_descriptions


# filter a list of star descriptions on the presence of a catalog
def catalog_filter(star: StarDescription, catalog_name):
    if star.get_catalog(catalog_name) is not None:
        return True
    return False


def interact():
    import code
    code.InteractiveConsole(locals=dict(globals(), **locals())).interact()


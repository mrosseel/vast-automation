from typing import List
from star_description import StarDescription
import do_charts_vast
from functools import partial
import toml
import os
from reading import trash_and_recreate_dir, create_dir
from shutil import copy
import glob
import logging
from pathlib import PurePath
import pytz
from datetime import datetime
import utils
from star_metadata import SiteData

UNKNOWN = "Unknown"


def run(post_name: str, selected_stars: List[StarDescription], len_vsx: int, len_candidates: int, resultdir: str):
    sitedir = f"{os.getcwd()}/site/vsx/"
    images_prefix = f"/images/{post_name}/"
    copy_files(post_name, resultdir, sitedir)
    result = get_header(post_name)
    result += get_starfile_preamble(images_prefix, len(selected_stars), len_vsx, len_candidates)
    sorted_stars = utils.sort_selected(selected_stars)
    part_block = partial(block, resultdir=resultdir, images_prefix=images_prefix)
    for star in sorted_stars:
        result += part_block(star)
    postdir = f"{sitedir}/content/posts/{post_name}/"
    create_dir(postdir)
    with open(f"{postdir}/{post_name}.md", 'w') as outfile:
        outfile.write(result)


def copy_files(post_name: str, resultdir: str, sitedir: str):
    imagesdir = f"{sitedir}static/images/{post_name}/"
    trash_and_recreate_dir(imagesdir)
    selected_phase = f'{resultdir}phase_selected/*.png'
    selected_phase_glob = glob.glob(selected_phase)
    vsx_phase = f'{resultdir}phase_vsx/*.png'
    vsx_phase_glob = glob.glob(vsx_phase)
    candidates_phase = f'{resultdir}phase_candidates/*.png'
    candidates_phase_glob = glob.glob(candidates_phase)
    aavso = f'{resultdir}aavso*/*.txt'
    aavso_glob = glob.glob(aavso)
    fieldcharts = f'{resultdir}fieldcharts/*.png'
    fieldcharts_glob = glob.glob(fieldcharts)
    lightcharts = f'{resultdir}light*/*.png'
    lightcharts_glob = glob.glob(lightcharts)
    logging.info(f"Copying {len(selected_phase_glob)} phase files from {selected_phase}...")
    for file in selected_phase_glob:
        copy(file, imagesdir)
    logging.info(f"Copying {len(vsx_phase_glob)} phase files from {vsx_phase}...")
    for file in vsx_phase_glob:
        copy(file, imagesdir)
    logging.info(f"Copying {len(candidates_phase_glob)} phase files from {candidates_phase}...")
    for file in candidates_phase_glob:
        copy(file, imagesdir)
    logging.info(f"Copying {len(lightcharts_glob)} light charts from {lightcharts}...")
    for file in lightcharts_glob:
        copy(file, imagesdir)
    logging.info(f"Copying {len(aavso_glob)} aavso files from {aavso}...")
    for file in aavso_glob:
        copy(file, imagesdir)
    logging.info(f"Copying {len(fieldcharts_glob)} field charts from {fieldcharts}...")
    for file in fieldcharts_glob:
        copy(file, imagesdir)
    copy(f"{resultdir}selected_radec.txt", imagesdir)
    logging.info(f"Copying selected_radec.txt...")
    logging.info(f"Copying done.")


def block(star: StarDescription, resultdir: str, images_prefix: str):
    try:
        is_vsx = star.has_metadata("VSX")
        is_selected = star.has_metadata("SELECTEDFILE")
        is_candidate = star.has_metadata("CANDIDATE")
        sitedata: SiteData = star.get_metadata("SITE")

        _, _, extradata, filename_no_ext = utils.get_star_or_catalog_name(star, suffix="")
        filename_no_ext_phase = filename_no_ext + "_phase"
        if is_vsx:
            txt_path = PurePath(resultdir, 'phase_vsx/txt', filename_no_ext + '.txt')
        elif is_candidate:
            txt_path = PurePath(resultdir, 'phase_candidates/txt', filename_no_ext + '.txt')
        elif is_selected:
            txt_path = PurePath(resultdir, 'phase_selected/txt', filename_no_ext + '.txt')

        try:
            parsed_toml = toml.load(txt_path)
        except FileNotFoundError:
            logging.error(f"Could not load txt file with phase information from {txt_path}")

        ucac4 = star.get_metadata("UCAC4")
        if ucac4 is None:
            ucac4_name = f"{star.coords}"
        else:
            ucac4_name = ucac4.catalog_id if not None else UNKNOWN
        nl = '\n'
        name = f"{nl.join(parsed_toml['our_name'])}" if 'our_name' in parsed_toml else f"OUR_NAME_{star.local_id}"
        period = f"{float(parsed_toml['period']):.5f}"
        phase_url = f"{images_prefix}{filename_no_ext_phase}.png"
        epoch = f"{parsed_toml['epoch']}" if 'epoch' in parsed_toml else UNKNOWN
        var_type = f"{parsed_toml['var_type']}" if 'var_type' in parsed_toml else UNKNOWN
        vsx_var_flag = f" ({parsed_toml['vsx_var_flag']})" if 'vsx_var_flag' in parsed_toml else ""
        vsx_separation = f"<li>VSX separation: +/- {parsed_toml['vsx_separation']*3600:.0f} arcsec</li>" \
            if 'vsx_separation' in parsed_toml else ""
        var_type_link = f"<a href='https://www.aavso.org/vsx/index.php?view=help.vartype&nolayout=1&abbrev=" \
                        f"{var_type}'>{var_type}</a>" if var_type != UNKNOWN else var_type
        mag_range = f"{parsed_toml['range']}"
        minmax = f"<li>minmax: {parsed_toml['minmax']}</li>" if 'minmax' in parsed_toml else ""
        vsx_link = f'<li><a target="_blank" rel="noopener noreferrer" ' \
                   f'href="https://www.aavso.org/vsx/index.php?view=detail.top&oid={extradata["OID"]}"' \
                   f'>VSX link</a></li>' if is_vsx else ""
        points_removed = f"<li>Outliers removed: {parsed_toml['points_removed']}</li>" \
            if parsed_toml['points_removed'] > 0 else ""
        result = f'''<div class="bb-l b--black-10 w-100">
        <div class="fl w-70 pa2 ba">
            <img class="special-img-class" src="{phase_url}" alt="{phase_url}"/>
        </div>
        <div class="fl w-30 pa2 ba">
            <ul>
            <li>{ucac4_name}</li>
            <li>name: {name}</li>
            <li>period (d): {period}</li>{minmax}
            <li>mag. range: {mag_range}</li>
            <li><a target="_blank" rel="noopener noreferrer" href="
            https://www.aavso.org/vsx/index.php?view=about.vartypessort">type</a>: {var_type_link}{vsx_var_flag}</li>{vsx_link}{vsx_separation}
            <li>epoch: {epoch}</li>
            <li>coords: {utils.get_hms_dms(star.coords)}</li>{points_removed}
            <li><a href="{images_prefix}vsx_and_star_{star.local_id:05}.png">finder chart</a></li>
            <li><a href="{images_prefix}{filename_no_ext}_ext.txt">observations</a></li>
            <li><a href="{images_prefix}{filename_no_ext}_light.png">light curve</a></li>
            <li><a href="{images_prefix}{filename_no_ext}_lightpa.png">light curve PA</a></li>
            </ul>
        </div>
    </div>
    '''
        return result
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        import traceback
        print(traceback.print_exc())
        logging.error(message)
        logging.error("File not found error in store and curve for star", star.path)
        return f'<div class="fl w-100 pa2 ba">Could not load {txt_path}</div>'


def get_header(title: str):
    return f'---\ntitle: "{title}"\ndate: {get_date_time()}\ndraft: false\n' \
           f'summary: "Batch of new variable stars, created at {datetime.now()}"\n---\n'


def get_date_time():
    return pytz.utc.localize(datetime.utcnow()).isoformat()


def get_starfile_preamble(images_prefix: str, len_selected: int, len_vsx: int, len_candidates: int):
    return f'<div class="bb-l b--black-10 w-100">' \
           f'<div class="fl w-100 pa2 ba">' \
           f'Images by Josch Hambsch, Data processing by Mike Rosseel, Josch Hambsch and ' \
           f'<a href="http://scan.sai.msu.ru/vast/">VaST</a></div>' \
           f'<a href="{images_prefix}selected_radec.txt">CSV file of all stars on this page</a><br/>' \
           f'<a href="{images_prefix}vsx_{len_vsx}_and_selected_{len_selected}.png">' \
           f'Finder chart with VSX and selected new stars</a><br/>' \
           f'Periods are derived using Lomb-Scargle (LS), Peranso (OWN) or from the VSX database (VSX)</div></div>\n'

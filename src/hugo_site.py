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
from star_metadata import StarFileData


def run(post_name: str, selected_stars: List[StarDescription], resultdir: str):
    sitedir = f"{os.getcwd()}/site/vsx/"
    copy_files(post_name, resultdir, sitedir)
    result = get_header(post_name)
    sorted_stars = sorted(selected_stars, key=lambda x: x.get_metadata('STARFILE').our_name)
    part_block = partial(block, resultdir=resultdir, post_name=post_name)
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
    aavso = f'{resultdir}aavso*/*.txt'
    aavso_glob = glob.glob(aavso)
    fieldcharts = f'{resultdir}fieldcharts/*.png'
    fieldcharts_glob = glob.glob(fieldcharts)
    lightcharts = f'{resultdir}light*/*.png'
    lightcharts_glob = glob.glob(lightcharts)
    logging.info(f"Copying {len(selected_phase_glob)} phase files from {selected_phase}...")
    for file in selected_phase_glob:
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
    logging.info(f"Copying done.")


def block(star: StarDescription, resultdir: str, post_name: str):
    try:
        vsx_name, separation, filename_no_ext = do_charts_vast.get_star_or_vsx_name(star, suffix="_phase")
        txt_path = PurePath(resultdir, 'phase_candidates/txt', filename_no_ext + '.txt')
        metadata: StarFileData = star.get_metadata("STARFILE")
        try:
            parsed_toml = toml.load(txt_path)
        except FileNotFoundError:
            # txt can be in candidates or selected.
            txt_path = PurePath(resultdir, 'phase_selected/txt', filename_no_ext + '.txt')
            parsed_toml = toml.load(txt_path)
        ucac4 = star.get_metadata("UCAC4")
        if ucac4 is None:
            ucac4_name = f"{star.coords}"
        else:
            ucac4_name = ucac4.catalog_id if not None else "unknown"
        period = f"{parsed_toml['period']:.5f}" if metadata.period is None \
            else f"{metadata.period:.5f} +/- {metadata.period_err:.5f}"
        images_prefix = f"/images/{post_name}/"
        phase_url = f"{images_prefix}{filename_no_ext}.png"
        minmax = f"<li>{metadata.minmax}</li>" if metadata.minmax is not None else ""
        epoch = f"<li>epoch: {metadata.epoch}</li>" if metadata.epoch is not None else ""
        var_type = metadata.var_type if metadata.var_type else "Unknown"
        result = f'''<div class="bb-l b--black-10 w-100">
        <div class="fl w-70 pa2 ba">
            <img class="special-img-class" src="{phase_url}" alt="{phase_url}"/>
        </div>
        <div class="fl w-30 pa2 ba">
            <ul>
            <li>{ucac4_name}</li>
            <li>name: {metadata.our_name}</li>
            <li>period (d): {period}</li>
            {minmax}
            <li>mag. range: {parsed_toml['range']}</li>
            <li>type: {var_type}</li>
            {epoch}
            <li>coords: {utils.get_hms_dms(star.coords)}</li>
            <li><a href="{images_prefix}vsx_and_star_{star.local_id:05}.png">finder chart</a></li>
            <li><a href="{images_prefix}{star.local_id:05}_ext.txt">observations</a></li>
            <li><a href="{images_prefix}{star.local_id:05}_light.png">light curve</a></li>
            </ul>
        </div>
    </div>
    '''
        return result
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        logging.error(message)
        logging.error("File not found error in store and curve for star", star.path)
        return f'<div class="fl w-100 pa2 ba">Could not load {txt_path}</div>'


def get_header(title: str):
    return f'---\ntitle: "{title}"\ndate: {pytz.utc.localize(datetime.utcnow()).isoformat()}\ndraft: false\n' \
           f'summary: "Batch of new variable stars"\n---\n'

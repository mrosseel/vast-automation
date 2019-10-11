from typing import List
from star_description import StarDescription
from functools import partial
import toml
import os
from reading import trash_and_recreate_dir, create_dir
from shutil import copy
import glob
import logging


def run(post_name: str, selected_stars: List[StarDescription], resultdir: str):
    sitedir = f"{os.getcwd()}/site/vsx/"
    copy_files(post_name, resultdir, sitedir)
    result = get_header(post_name)
    part_block = partial(block, resultdir=resultdir, post_name=post_name)
    for star in selected_stars:
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
    aavso = f'{resultdir}aavso/*.txt'
    aavso_glob = glob.glob(aavso)
    fieldcharts = f'{resultdir}fieldcharts/*.png'
    fieldcharts_glob = glob.glob(fieldcharts)
    logging.info(f"Copying {len(selected_phase_glob)} phase files from {selected_phase}...")
    for file in selected_phase_glob:
        copy(file, imagesdir)
    logging.info(f"Copying {len(aavso_glob)} aavso files from {aavso}...")
    for file in aavso_glob:
        copy(file, imagesdir)
    logging.info(f"Copying {len(fieldcharts_glob)} field charts from {fieldcharts}...")
    for file in fieldcharts_glob:
        copy(file, imagesdir)


def block(star: StarDescription, resultdir: str, post_name: str):
    try:
        star_match, separation = star.get_match_string("VSX")
        phase_file = f"{star_match}_phase.png" if star_match is not None else f"{star.local_id:05}_phase.png"
        parsed_toml = toml.load(f"{resultdir}/phase_candidates/txt/{phase_file[:-4]}.txt")
        ucac4 = star.get_catalog("UCAC4", strict=True)
        ucac4_name = ucac4 if not None else "unknown"
        images_prefix = f"/images/{post_name}/"
        phase_url = f"{images_prefix}{phase_file}"
        result = f'''<div class="bb-l b--black-10 w-100">
        <div class="fl w-70 pa2 ba">
            <img class="special-img-class" src="{phase_url}" alt="{phase_url}"/>
        </div>
        <div class="fl w-30 pa2 ba">
            <ul>
            <li>UCAC4: {ucac4_name}</li>
            <li>period (d): {parsed_toml['period']:.5f}</li>
            <li>mag. range: {parsed_toml['range']}</li>
            <li>type: </li>
            <li>coords: {parsed_toml['coords'][0]:.5f} {parsed_toml['coords'][1]:.5f}</li>
            <li><a href="{images_prefix}vsx_and_star_{star.local_id:05}.png">finder chart</a></li>
            <li><a href="{images_prefix}{star.local_id:05}_ext.txt">aavso observations</a></li>
            
            </ul>
        </div>
    </div>
    '''
        return result
    except:
        return f'<div class="fl w-100 pa2 ba">Could not load {phase_file}</div>'


def get_header(title: str):
    return f'---\ntitle: "{title}"\ndate: 2019-10-09T11:10:29+02:00\ndraft: false\n' \
           f'summary: "Batch of new variable stars"\n---\n'

from typing import List
from star_description import StarDescription
from pathlib import Path
from functools import partial
import toml
import os
from reading import trash_and_recreate_dir, create_dir
from shutil import copy
import glob
import logging
import pytz
from datetime import datetime
import utils
from star_metadata import SiteData

UNKNOWN = "Unknown"


def run(
    post_name: str,
    selected_stars: List[StarDescription],
    len_vsx: int,
    reference_frame: Path,
    resultdir: str,
):
    sitedir = f"{os.getcwd()}/site/vsx/"
    images_prefix = f"/images/{post_name}/"
    # copy_files(post_name, resultdir, sitedir)
    staticimagesdir = f"{sitedir}static/{images_prefix}"
    selective_copy_files(selected_stars, staticimagesdir, resultdir)
    result = get_header(post_name, reference_frame.name)
    result += get_starfile_preamble(
        images_prefix,
        len([x for x in selected_stars if not x.has_metadata("VSX")]),
        len_vsx,
        get_optional_preamble(images_prefix, staticimagesdir),
    )
    sorted_stars = utils.sort_selected(selected_stars)
    part_block = partial(block, resultdir=resultdir, images_prefix=images_prefix)
    for star in sorted_stars:
        result += part_block(star)
    postdir = f"{sitedir}/content/posts/{post_name}/"
    create_dir(postdir)
    with open(f"{postdir}/{post_name}.md", "w") as outfile:
        outfile.write(result)


def get_optional_preamble(images_prefix: str, destdir: str):
    ap_air = "aperture_vs_airmass.png"
    ap_jd = "aperture_vs_jd.png"
    ap_air_img = Path(destdir, ap_air)
    ap_jd_img = Path(destdir, ap_jd)
    optionalpreamble = ""
    if ap_air_img.is_file() and ap_jd_img.is_file():
        optionalpreamble = f'<a href="{images_prefix}/{ap_air}" alt="Plot of aperture vs airmass">Aperture vs Airmass</a> and '
        optionalpreamble += f'<a href="{images_prefix}/{ap_jd}" alt="Plot of aperture vs Julian Day">Aperture vs JD</a><br/>'
    return optionalpreamble + "<br/>"


def selective_copy_files(stars: List[StarDescription], destdir: str, resultdir: str):
    trash_and_recreate_dir(destdir)
    for astar in stars:
        copy(astar.result["phase"], destdir) if "phase" in astar.result else None
        copy(
            astar.result["compstars"], destdir
        ) if "compstars" in astar.result else None
        copy(astar.result["compA"], destdir) if "compA" in astar.result else None
        copy(astar.result["compB"], destdir) if "compB" in astar.result else None
        copy(astar.result["light"], destdir) if "light" in astar.result else None
        copy(astar.result["lightpa"], destdir) if "lightpa" in astar.result else None
        copy(astar.result["lightcont"], destdir) if "lightcont" in astar.result else None
        copy(astar.result["lightmain"], destdir) if "lightmain" in astar.result else None
        copy(astar.result["aavso"], destdir) if "aavso" in astar.result else None
    fieldcharts = f"{resultdir}fieldcharts/*.png"
    fieldcharts_glob = glob.glob(fieldcharts)
    logging.info(f"Copying {len(fieldcharts_glob)} field charts from {fieldcharts}...")
    for file in fieldcharts_glob:
        copy(file, destdir)
    copy(f"{resultdir}selected_radec.txt", destdir)
    copy(f"{resultdir}selected_localid.txt", destdir)
    copy(f"{resultdir}aavso_vsx.txt", destdir)
    logging.info(f"Copying done.")

def get_from_toml(key, parsed_toml, default=None):
    if key in parsed_toml:
        return parsed_toml[key]
    else:
        return default

def block(star: StarDescription, resultdir: str, images_prefix: str):
    try:
        is_vsx = star.has_metadata("VSX")
        starui: utils.StarUI = utils.get_star_or_catalog_name(star, suffix="_phase")
        txt_path = Path(
            Path(star.result["phase"]).parent,
            "txt",
            starui.filename_no_suff_no_ext + ".txt",
        )
        try:
            parsed_toml = toml.load(txt_path)
        except FileNotFoundError:
            logging.error(
                f"Could not load txt file with phase information from {txt_path}"
            )

        ucac4 = star.get_metadata("UCAC4")
        if ucac4 is None:
            ucac4_name = f"no UCAC4 match !!!"
            ucac4_mag = f"no mag"
            ucac4_coords = f""

        else:
            ucac4_name = ucac4.catalog_id if not None else UNKNOWN
            ucac4_mag = (
                f'<strong style="color: red;">{ucac4.vmag}</strong>'
                if abs(ucac4.vmag - parsed_toml["vmag"]) > 0.2 and ucac4.vmag != 20.0
                else ucac4.vmag
            )
            ucac4_coords = (
                f"<li>coords: {utils.get_hms_dms_sober(ucac4.coords)} (ucac4)</li>"
            )
        nl = "\n"
        name = (
            f"\n{parsed_toml['our_name']}" if "our_name" in parsed_toml else f"OUR_NAME_{star.local_id}"
        )
        period = f"{float(parsed_toml['period']):.5f}" if 'period' in parsed_toml else 'None'
        var_type_raw = get_from_toml('var_type', parsed_toml, UNKNOWN)
        var_type = f"{var_type_raw}"
        phase_url = f"{images_prefix}{starui.filename_no_ext}.png"
        if utils.is_var_type_periodic(var_type):
            phase_url = f"{images_prefix}{starui.filename_no_suff_no_ext}_lightmain.png"
        epoch = f"{parsed_toml['epoch']}" if 'epoch' in parsed_toml else UNKNOWN
        vsx_var_flag = f" ({parsed_toml['vsx_var_flag']})" if 'vsx_var_flag' in parsed_toml else ""
        tomlseparation = parsed_toml['separation'] if 'separation' in parsed_toml else None
        ucacseparation = star.coords.separation(star.get_metadata("UCAC4").coords).degree if star.has_metadata(
            "UCAC4") else None
        realseparation = ucacseparation if ucacseparation is not None else tomlseparation if tomlseparation is not None else None
        separation = f"<li>separation: +/- {realseparation * 3600:.0f} arcsec</li>" if realseparation is not None else ""
        var_type_link = f"<a href='https://www.aavso.org/vsx/index.php?view=help.vartype&nolayout=1&abbrev=" \
                        f"{var_type}'>{var_type}</a>" if var_type != UNKNOWN else var_type
        mag_range = f"{parsed_toml['range']}"
        minmax = (
            f"<li>minmax: {parsed_toml['minmax']}</li>"
            if "minmax" in parsed_toml
            else ""
        )
        vsx_link = (
            f'<li><a target="_blank" rel="noopener noreferrer" '
            f'href="https://www.aavso.org/vsx/index.php?view=detail.top&oid={starui.extradata["OID"]}"'
            f">VSX link</a></li>"
            if is_vsx
            else ""
        )
        points_removed = (
            f"<li>Outliers removed: {parsed_toml['points_removed']}</li>"
            if parsed_toml["points_removed"] > 0
            else ""
        )
        optional_compstars = (
            f'<a href="{images_prefix}{starui.filename_no_suff_no_ext}_compstarsA.png" '
            f'alt="Plot of all comparison stars used to measure this star">C</a>, '
            f'<a href="{images_prefix}{starui.filename_no_suff_no_ext}_compstarsB.png" '
            f'alt="Plot of all comparison stars used to measure this star + the star itself">C+V</a>, '
            if "compA" in star.result
            else ""
        )
        optional_stats = (
            f'<li>stats: <a href="{images_prefix}{starui.filename_no_suff_no_ext}_merr_vs_jd.png"  '
            f'alt="Plot of magnitude error vs JD">merr_vs_jd</a></li>'
            if "merr_vs_jd" in star.result
            else ""
        )
        result = f"""<div class="bb-l b--black-10 w-100">
        <div class="fl w-70 pa2 ba">
            <img class="special-img-class" src="{phase_url}" alt="{phase_url}"/>
        </div>
        <div class="fl w-30 pa2 ba">
            <ul>
            <li>{ucac4_name} (mag:{ucac4_mag})</li>
            <li>name: {name}</li>{ucac4_coords}
            <li>coords: {utils.get_hms_dms_sober(star.coords)} (ours)</li>{separation}{points_removed}
            <li>period (d): {period}</li>{minmax}
            <li>mag. range: {mag_range}</li>
            <li><a target="_blank" rel="noopener noreferrer" href="https://www.aavso.org/vsx/index.php?view=about.vartypessort">type</a>: {var_type_link}{vsx_var_flag}</li>
            {vsx_link}<li>epoch: {epoch}</li>
            <li><a href="{images_prefix}vsx_and_star_{starui.filename_no_suff_no_ext}.png">finder chart</a></li>
            <li><a href="{images_prefix}{starui.filename_no_suff_no_ext}_ext.txt">observations</a></li>
            <li>light curve: <a href="{images_prefix}{starui.filename_no_suff_no_ext}_light.png">Normal</a>,
            <a href="{images_prefix}{starui.filename_no_suff_no_ext}_lightpa.png">PA</a>,
            <a href="{images_prefix}{starui.filename_no_suff_no_ext}_lightcont.png">Continuous</a></li>
            <li>comparison stars: {optional_compstars}<a href="{images_prefix}{starui.filename_no_suff_no_ext}_comps.txt">list</a></li>{optional_stats}
            </ul>
        </div>
    </div>
    """
        return result
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        import traceback

        print(traceback.print_exc())
        logging.error(message)
        logging.error("File not found error in store and curve for star" + star.path)
        return f'<div class="fl w-100 pa2 ba">Could not load {txt_path}</div>'


def get_header(title: str, ref_frame: str):
    return (
        f'---\ntitle: "{title}"\ndate: {get_date_time()}\ndraft: false\n'
        f'summary: "Using ref frame `{ref_frame}` - created at {datetime.now()}"\n---\n'
    )


def get_date_time():
    return pytz.utc.localize(datetime.utcnow()).isoformat()


def get_starfile_preamble(
    images_prefix: str, len_selected: int, len_vsx: int, optionalpreamble: str
):
    return (
        f'<div class="bb-l b--black-10 w-100">'
        f'<div class="fl w-100 pa2 ba">'
        f"Images by Josch Hambsch, Data processing by Mike Rosseel, Josch Hambsch and "
        f'<a href="http://scan.sai.msu.ru/vast/">VaST</a></div>'
        f'<a href="{images_prefix}aavso_vsx.txt">Output for AAVSO VSX</a><br/>'
        f'<a href="{images_prefix}selected_radec.txt">Internal use: csv with ra/dec of selected stars</a><br/>'
        f'<a href="{images_prefix}selected_localid.txt">Internal use: csv with vast ids of selected stars</a><br/>'
        f'<a href="{images_prefix}vsx_{len_vsx}_and_selected_{len_selected}.png">'
        f"Finder chart with {len_vsx} VSX and {len_selected} new stars</a><br/>{optionalpreamble}"
        f"Periods are derived using Lomb-Scargle (LS), Peranso (OWN) or from the VSX database (VSX)</div></div>\n"
    )

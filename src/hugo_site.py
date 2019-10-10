from typing import List
from star_description import StarDescription
from functools import partial
import toml

header = """---
title: "WW CrA new variables"
date: 2019-10-09T11:10:29+02:00
draft: true
---
"""


def run(selected_stars: List[StarDescription], resultdir: str):
    result = header
    part_block = partial(block, resultdir=resultdir)
    for star in selected_stars:
        result += part_block(star)
    with open(f"{resultdir}output.md", 'w') as outfile:
        outfile.write(result)


def block(star: StarDescription, resultdir: str):
    try:
        star_match, separation = star.get_match_string("VSX")
        phase_file = f"{star_match}_phase.png" if star_match is not None else f"{star.local_id:05}_phase.png"
        parsed_toml = toml.load(f"{resultdir}/phase_candidates/txt/{phase_file[:-4]}.txt")
        result = f'''<div class="bb-l b--black-10 w-100">
        <div class="fl w-70 pa2 ba">
            <img class="special-img-class" src="/munipack-automation/images/{phase_file}" alt="/munipack-automation/images/{phase_file}"/>
        </div>
        <div class="fl w-30 pa2 ba">
            <ul>
            <li>period: {parsed_toml['period']}</li>
            <li>range: {parsed_toml['range']}</li>
            <li>coords: {parsed_toml['coords'][0]} {parsed_toml['coords'][1]}</li>
            <li><a href="/munipack-automation/images/vsx_and_star_{star.local_id:05}.png">finder chart</a></li>
            <li><a href="/munipack-automation/images/vsx_and_star_{star.local_id:05}.png">finder chart</a></li>
            </ul>
        </div>
    </div>
    '''
        return result
    except:
        print("could not load", phase_file)
        return ""

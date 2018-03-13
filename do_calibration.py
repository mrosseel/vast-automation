import os
import init
import reading
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
import requests

def calibrate():
    w = WCS(init.reference_header)
    #xpos, ypos = w.all_world2pix(init.ra_deg, init.dec_deg, 0, ra_dec_order=True)
    #print(xpos, ypos)
    #print(w)
    return w

def find_reference_in_files(the_path):
    the_dir = os.listdir(the_path)
    the_dir.sort()
    reference_frame_index = the_dir.index(init.reference_frame)
    return reference_frame_index

def find_target_star(target_ra_deg, target_dec_deg, nr_results):
    target = SkyCoord(target_ra_deg, target_dec_deg, unit='deg')
    result_dict = reading.read_world_positions(init.worldposdir)
    distances_dict = {}
    for key in result_dict:
        distances_dict[key] = target.separation(SkyCoord(result_dict[key][0], result_dict[key][1], unit='deg')).degree
    df = pd.DataFrame(list(distances_dict.items()), columns=['star_nr', 'deg_separation'])
    df.sort_values(by='deg_separation', inplace=True)
    return df[:nr_results]

def find_target_stars(max_deg_separation):
    target = SkyCoord(target_ra_deg, target_dec_deg, unit='deg')
    result_dict = reading.read_world_positions(init.worldposdir)
    distances_dict = {}
    for key in result_dict:
        distances_dict[key] = target.separation(SkyCoord(result_dict[key][0], result_dict[key][1], unit='deg')).degree
    df = pd.DataFrame(list(distances_dict.items()), columns=['filename', 'deg_separation'])
    df.sort_values(by='deg_separation', inplace=True)
    return df[:nr_results]

def getVSX():
    payload = {'key1': 'value1', 'key2': 'value2'}
    r = requests.get('https://api.github.com/user', params=payload)

# returns [star_id, label, probability, flag, SkyCoord]
def getCandidates():
    df = pd.DataFrame.from_csv(init.basedir+'upsilon_output.txt')
    df.sort_values(by='probability', ascending=False)
    df=df[df['label'] != 'NonVar']
    df=df[df["probability"] > 0.50]
    df=df[df["flag"] != 1]
    positions=reading.read_world_positions(init.worldposdir)
    result = []
    for index, row in df.iterrows():
        #print(index, row)
        result.append([index, row['label'], row['probability'], row['flag'], SkyCoord(positions[int(index)][0], positions[int(index)][1], unit='deg')])
    return result

# returns [index, skycoord, type]
def getVSX():
    result = []
    df = pd.DataFrame.from_csv(init.basedir+'SearchResults.csv')
    #print(df.head())
    for index, row in df.iterrows():
        skycoord = SkyCoord(row['Coords'], unit=(u.hourangle, u.deg))
        result.append([index, skycoord, row['Type']])
    return result

# returns {'star_id': [label, probability, flag, SkyCoord, match_name, match_skycoord, match_type, separation_deg]}
def findNames():
    vsx = getVSX()
    candidates = getCandidates()
    result = {}
    for candidate in candidates:
        best_sep = 360
        best_var = None
        for variable in vsx:
            sep_deg = candidate[4].separation(variable[1]).degree
            if(sep_deg < best_sep):
                best_sep = sep_deg
                best_var = variable
        #candidate.append([best_var[0],best_var[1],best_var[2],best_sep])
        result[candidate[0]] = [candidate[1],candidate[2],candidate[3],candidate[4],best_var[0],best_var[1],best_var[2],best_sep]
    return result

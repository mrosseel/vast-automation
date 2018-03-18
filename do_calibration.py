import os
import init
import reading
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import numpy as np
import requests
from gatspy.periodic import LombScargleFast
import matplotlib as mp
mp.use('Agg') # needs no X server
import matplotlib.pyplot as plt

def calibrate():
    w = WCS(init.reference_header)
    #xpos, ypos = w.all_world2pix(init.ra_deg, init.dec_deg, 0, ra_dec_order=True)
    #print(xpos, ypos)
    #print(w)
    return w

def find_reference_frame_index():
    the_dir = os.listdir(init.fitsdir)
    the_dir.sort()
    reference_frame_index = the_dir.index(init.reference_frame)
    return reference_frame_index

# returns 'path + phot????.pht', the photometry file matched with the reference frame
def find_reference_photometry(reference_frame_index):
    the_dir = os.listdir(init.photometrydir)
    the_dir.sort()
    return init.photometrydir + the_dir[reference_frame_index]

def find_reference_matched(reference_frame_index):
    the_dir = os.listdir(init.matchedphotometrydir)
    the_dir.sort()
    return init.matchedphotometrydir + the_dir[reference_frame_index]

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

# returns [star_id, label, probability, flag, SkyCoord]
def getCandidates(threshold_prob=0.5):
    df = pd.DataFrame.from_csv(init.basedir+'upsilon_output.txt')
    df.sort_values(by='probability', ascending=False)
    df=df[df['label'] != 'NonVar']
    df=df[df["probability"] > threshold_prob]
    df=df[df["flag"] != 1]
    positions=reading.read_world_positions(init.worldposdir)
    result = []
    for index, row in df.iterrows():
        #print(index, row)
        result.append([index, row['label'], row['probability'], row['flag'], SkyCoord(positions[int(index)][0], positions[int(index)][1], unit='deg')])
    return result

# returns [index, skycoord, type]
def getVSX(the_file):
    result = []
    df = pd.DataFrame.from_csv(the_file)
    #print(df.head())
    for index, row in df.iterrows():
        skycoord = SkyCoord(row['Coords'], unit=(u.hourangle, u.deg))
        result.append([index, skycoord, row['Type']])
    return result

# returns {'star_id': [label, probability, flag, SkyCoord, match_name, match_skycoord, match_type, separation_deg]}
def find_vsx_for_upsilon_candidates(threshold_prob_candidates=0.5, max_separation=0.01):
    vsx = getVSX(init.basedir+'SearchResults.csv')
    candidates = getCandidates(threshold_prob_candidates)
    print("Got", len(candidates), "candidates and", len(vsx), "stars to check against.")
    result = {}
    for candidate in candidates:
        best_sep_deg = 360
        best_sep_string = ""
        best_var = None
        for variable in vsx:
            sep = candidate[4].separation(variable[1])
            if(sep.degree < max_separation and sep.degree < best_sep_deg):
                best_sep_deg = sep.degree
                best_sep_string = sep.to_string()
                best_var = variable
        if best_var == None:
            result[candidate[0]] = [candidate[1],candidate[2],candidate[3],candidate[4],'','','',best_sep_deg, best_sep_string]
        else:
            result[candidate[0]] = [candidate[1],candidate[2],candidate[3],candidate[4],best_var[0],best_var[1],best_var[2],best_sep_deg, best_sep_string]
    return result

# Takes in a list of known variables and maps them to the munipack-generated star numbers
# usage:
# vsx = getVSX(init.basedir+'SearchResults.csv')
# detections = reading.read_world_positions(init.worldposdir)
# returns { 'name of VSX variable': [VSX_var_SkyCoord, best_starfit, best_separation] }
def find_star_for_known_vsx(vsx, detections, max_separation=0.01):
    ra2 = np.array([])
    dec2 = np.array([])
    print("Creating astropy Catalog...")
    # rewrite ra/dec to skycoord objects
    for key in detections:
        ra2 = np.append(ra2, [detections[key][0]])
        dec2 = np.append(dec2, [detections[key][1]])
    print(ra2.shape, dec2.shape)
    catalog = SkyCoord(ra=ra2, dec=dec2, unit='deg')
    result = {}
    print("Searching best matches with max separation", max_separation, "...")
    for variable in vsx:
        print("Searching for", variable[0])
        idx, d2d, d3d = variable[1].match_to_catalog_sky(catalog)
        if d2d.degree < max_separation:
            result[variable[0]] = [variable[1], idx+1, d2d.degree]
            print("Found result for ", variable[0], ":", result[variable[0]])
    print("Found", len(result), "matches")
    return result

def calculate_phase_diagram(star):
    print("Calculating phase diagram for", star)
    curve = reading.read_lightcurve(star)
    if curve is None:
        return
    t_np = curve['JD'].as_matrix()
    y_np = curve['V-C'].as_matrix()
    dy_np = curve['s1'].as_matrix()
    ls = LombScargleFast()
    period_max = np.max(t_np)-np.min(t_np)
    ls.optimizer.period_range = (0.01,period_max)
    ls.fit(t_np,y_np)
    period = ls.best_period
    print("Best period: " + str(period) + " days")
    fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    plt.xlabel("Phase")
    plt.ylabel("Diff mag")
    plt.title("Lomb-Scargle-Periodogram for star "+str(star)+" period: " + str(period))

    # plotting + calculation of 'double' phase diagram from -1 to 1
    phased_t = np.fmod(t_np/period,1)
    minus_one = lambda t: t - 1
    minus_oner = np.vectorize(minus_one)
    phased_t2 = minus_oner(phased_t)
    phased_lc = y_np[:]
    phased_t_final = np.append(phased_t2, phased_t)
    phased_lc_final = np.append(phased_lc, phased_lc)
    phased_err = np.append(dy_np, dy_np)
    plt.errorbar(phased_t_final,phased_lc_final,yerr=phased_err,linestyle='none',marker='o')
    fig.savefig(init.phasedir+'phase'+str(star).zfill(5))
    plt.close(fig)

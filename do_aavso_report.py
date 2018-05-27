import aavso
import reading
from star_description import get_match_string
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
#supress the warning about vector transforms so as not to clutter the doc build log
import warnings
warnings.filterwarnings('ignore',module='astropy.coordinates.baseframe')
import init

def calculate_airmass(coord, location, jd):
    time = Time(jd, format='jd')
    altazs = coord.transform_to(AltAz(obstime=time, location=location))
    return altazs.secz

def report(star_description, comparison_star):
    #         row.append(observation_data['name'])
    #         row.append(observation_data['date'])
    #         row.append(observation_data['magnitude'])
    #         row.append(observation_data['magnitude_error'])
    #         row.append(observation_data['filter'])
    #         row.append(observation_data['transformed'])
    #         row.append(observation_data['magnitude_type'])
    #         row.append(observation_data.get('comparison_name', 'na'))
    #         row.append(observation_data.get('comparison_magnitude', 'na'))
    #         row.append(observation_data.get('check_name', 'na'))
    #         row.append(observation_data.get('check_magnitude', 'na'))
    #         row.append(observation_data.get('airmass', 'na'))
    #         row.append(observation_data.get('group', 'na'))
    #         row.append(observation_data.get('chart', 'na'))
    #         row.append(observation_data.get('notes', 'na'))
    curve = reading.read_lightcurve(star_description.local_id, filter=True, preprocess=False)
    star_match, separation = get_match_string(star_description)
    print(comparison_star)
    comparison_star_vmag = comparison_star[0].vmag
    title = str(star_description.local_id if star_description.aavso_id is None else star_description.aavso_id)
    earth_location = EarthLocation(lat=init.sitelat, lon=init.sitelong, height=390*u.m) # TODO WRONG !!!!!!!!!!!!

    with open(init.aavso_reports + title+'_extended.txt', 'w') as fp:
        writer = aavso.ExtendedFormatWriter(fp, 'RMH', software='munipack-automation', obstype='CCD')
        for index, row in curve.iterrows():
            #print(row, type(row))
            writer.writerow({
                'name': star_match if not star_match == '' else 'Unknown',
                'date': row['JD'],
                'magnitude': row['V-C'] + comparison_star_vmag,
                'magnitude_error': row['s1'],
                'filter': 'V',
                'transformed': 'YES',
                'magnitude_type': 'STD',
                'comparison_name': 'ENSEMBLE',
                'comparison_magnitude': 'na',
                'check_name': comparison_star[0].aavso_id if not comparison_star[0].aavso_id is None else '???',
                'check_magnitude': comparison_star_vmag,
                'airmass': calculate_airmass(star_description.coords, earth_location, row['JD']),
                'group': 'na',
                'chart': 'na',
                'notes': 'na'
            })

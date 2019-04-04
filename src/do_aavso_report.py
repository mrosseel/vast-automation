import aavso
import reading
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
#supress the warning about vector transforms so as not to clutter the doc build log
import warnings
warnings.filterwarnings('ignore',module='astropy.coordinates.baseframe')
from init_loader import init, settings
import tqdm

def calculate_airmass(coord, location, jd):
    time = Time(jd, format='jd')
    altazs = coord.transform_to(AltAz(obstime=time, location=location))
    return altazs.secz

def report(target_dir, star_description, comparison_star):
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
    star_match_ucac4, separation = star_description.get_match_string("UCAC4")
    star_match_vsx, separation = star_description.get_match_string("VSX", strict=False)
    comp_ucac4 = comparison_star.get_match_string("UCAC4", strict=True)
    var_display_name = star_match_ucac4 if star_match_vsx == None else star_match_vsx
    check_display_name = comparison_star.aavso_id if not comparison_star.aavso_id is None else comp_ucac4[0]

    print(" Star match:{}, comparison_star:{}".format(var_display_name, comparison_star))
    comparison_star_vmag = comparison_star.vmag
    title = str(star_description.local_id if star_description.aavso_id is None else star_description.aavso_id)
    earth_location = EarthLocation(lat=init.sitelat, lon=init.sitelong, height=init.sitealt*u.m)
    print("Starting aavso report with star:{}".format(star_description))

    with open(target_dir + title+'_extended.txt', 'w') as fp:
        writer = aavso.ExtendedFormatWriter(fp, 'RMH', software='munipack-automation', obstype='CCD')
        for _, row in tqdm.tqdm(curve.iterrows(), total=len(curve), unit="observations"):
            #print(row, type(row))
            writer.writerow({
                'name': var_display_name,
                'date': row['JD'],
                'magnitude': row['V-C'] + comparison_star_vmag,
                'magnitude_error': row['s1'],
                'filter': 'V',
                'transformed': 'YES',
                'magnitude_type': 'STD',
                'comparison_name': 'ENSEMBLE',
                'comparison_magnitude': 'na',
                'check_name': check_display_name,
                'check_magnitude': comparison_star_vmag,
                'airmass': calculate_airmass(star_description.coords, earth_location, row['JD']),
                'group': 'na',
                'chart': 'na',
                'notes': 'na'
            })

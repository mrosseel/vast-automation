import numpy as np
import multiprocessing

#aperture_range = np.arange(0.2, 10, 0.2)
aperture_range = np.arange(9, 10, 0.5)

all = range(1, 10000)
thousand = range(1, 1000)

#### SETTINGS for DO_MUNIWIN.PY ####

datadir = 'current/'
star_list = all #[52]
nr_threads = multiprocessing.cpu_count()*2

do_convert_fits=0
do_photometry=0
do_match=0
do_munifind=1
do_lightcurve=0
do_lightcurve_resume=0
do_pos=0
do_pos_resume=1
do_ml=0
do_lightcurve_plot=0
do_phase_diagram=0
do_field_charts=1
do_reporting=0

### CALIBRATION ###
reference_file = 'new-image.fits'
ra_deg = 271.4032917
dec_deg = -43.8326111
sitelat = '-22 57 10'
sitelong = '-68 10 49'
sitealt= 350 # TODO WRONG !!!! Used in AAVSO reporting for atmospheric extinction
### CALIBRATION


#### DO NOT TOUCH ####
# standard directories
codedir = './'
basedir = codedir + datadir
memdir = '/media/vars/'
reference_header = basedir + reference_file
# input + processing
fitsdir = basedir + "fits/"
convfitsdir = basedir + "converted_fits/"
photometrydir = basedir + "photometry/"
matchedphotometrydir = basedir + "matched_photometry/"
#matchedphotometrydir = memdir + "matched_photometry/"
aperturedir = basedir + "aperture/"
vsxcatalogdir = codedir + "vsx_catalog.bin"

# results
resultdir = basedir + "results/"
lightcurvedir = basedir + "lightcurves/"
posdir = resultdir + "positions/"
worldposdir = resultdir + "world_positions/"
chartsdir = resultdir + "charts/"
phasedir = resultdir + "phasediagrams/"
aavsoreportsdir = resultdir + "aavso_reports/"
fieldchartsdirs = resultdir + "field_charts/"

#### DO NOT TOUCH ####


################## CUSTOM star selections, can be ignored

selection = [58, 137, 138, 886]
selection2 = (1, 2, 3, 73, 138, 143, 264, 2675, 1045, 847, 1193)
custom_muniwin = range(5776, 10000)
# custom_charts = range(1,5776)
custom_charts = range(1, 2164)
wwcra_candidates = [6314,
                    6291,
                    5978,
                    5745,
                    5557,
                    5430,
                    4962,
                    3368,
                    2659,
                    2535,
                    2058,
                    1870,
                    1858,
                    1636,
                    1619,
                    1324,
                    464,
                    454,
                    441,
                    440,
                    435,
                    427,
                    420,
                    378,
                    363,
                    313,
                    290,
                    277,
                    227,
                    225,
                    205,
                    160,
                    159,
                    141,
                    132,
                    117,
                    111,
                    107,
                    1
                    ]
wwcra_certain_candidates = [6314,
                            5978,
                            5557,
                            5430,
                            4962,
                            3368,
                            2535,
                            2058,
                            1870,
                            1636,
                            464,
                            441,
                            427,
                            363,
                            290,
                            277,
                            227,
                            205,
                            141]

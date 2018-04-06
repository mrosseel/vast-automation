#### DO NOT TOUCH ####
import multiprocessing

codedir = './'

aperture = 2
all = range(1, 10000)
thou = range(1, 1000)
selection = [58, 137, 138, 886]
selection2 = (1, 2, 3, 73, 138, 143, 264, 2675, 1045, 847, 1193)
custom_muniwin = range(5776, 10000)
# custom_charts = range(1,5776)
custom_charts = range(1, 2164)
all_star_list = range(1, 138)
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

#### DO NOT TOUCH ####

#### SETTINGS for DO_MUNIWIN.PY ####

basedir = codedir + 'inputfiles/RTGru/'
#basedir = codedir + 'inputfiles/WWCrA_allflat/'
#basedir =  codedir + 'inputfiles/WWCrA2015/'
#basedir =  codedir + 'inputfiles/WWCrA_bigger/'
star_list = all
nr_threads = multiprocessing.cpu_count() 

do_convert_fits=False
do_photometry=False
do_match=False
do_munifind=False
do_lightcurve=False
do_lightcurve_resume=False
do_pos=False
do_pos_resume=False
do_calibrate=False
do_ml=False
do_charting=True
do_phase_diagram=True
do_reporting=False

#### DO NOT TOUCH ####
# standard directories
fitsdir = basedir + "fits/"
convfitsdir = basedir + "converted_fits/"
photometrydir = basedir + "photometry/"
matchedphotometrydir = basedir + "matched_photometry/"
resultdir = basedir + "results/"
lightcurvedir = resultdir + "lightcurves/"
posdir = resultdir + "positions/"
worldposdir = resultdir + "world_positions/"
chartsdir = resultdir + "charts/"
phasedir = resultdir + "phasediagrams/"
vsx_catalog_path = codedir + "vsx_catalog.bin"
#### DO NOT TOUCH ####

### CALIBRATION ###
reference_dir = convfitsdir
reference_frame = 'kout000546.fts'
reference_header = basedir + 'new-image.fits'
ra_deg = 271.4032917
dec_deg = -43.8326111
sitelat = '-22 57 10'
sitelong = '-68 10 49'
### CALIBRATION

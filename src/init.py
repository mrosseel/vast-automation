from multiprocessing import cpu_count

all = range(1, 10001)
thousand = range(1, 1000)

#### SETTINGS for DO_MUNIWIN.PY ####

datadir = 'current/'
aperture_find_percentage=0.05
free_memory_GB=10
star_list = all #[52]
nr_threads = cpu_count()

do_conf_phot=1
do_conf_match=0
do_convert_fits=0
do_photometry=0
photometry_wildcard=None  # init.convfitsdir+"kout0000??.fts"
do_match=1
do_aperture_search=1
do_lightcurve=1
do_lightcurve_resume=0
do_pos=1
do_pos_resume=0
do_ml=0
do_lightcurve_plot=1
do_phase_diagram=1
do_field_charts=1
do_reporting=0

### CALIBRATION ###
reference_file = 'new-image.fits'
ra_deg = 271.4032917
dec_deg = -43.8326111
sitelat = '-22 57 10'
sitelong = '-68 10 49'
sitealt= 2500
### CALIBRATION


#### DO NOT TOUCH ####
# standard directories
codedir = './'
basedir = codedir + datadir
reference_header = basedir + reference_file
# input + processing
fitsdir = basedir + "fits/"
convfitsdir = basedir + "converted_fits/"
photometrydir = basedir + "photometry/"
matchedphotometrydir = basedir + "matched_photometry/"
aperturedir = basedir + "aperture/"
vsxcatalogdir = codedir + "vsx_catalog.bin"

# results
resultdir = basedir + "results/"
lightcurvedir = resultdir + "lightcurves/"
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

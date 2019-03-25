from multiprocessing import cpu_count

all = range(1, 10001)
thousand = range(1, 1000)

#### SETTINGS for DO_MUNIWIN.PY ####

datadir = 'current/'
aperture_find_percentage=0.05
free_memory_GB=10
star_list = all
nr_threads = cpu_count()

# Comparison stars to use:
comparison_stars = ['UCAC4 231-154752', 'UCAC4 232-147677']
# use None if you want automatic aperture
aperture = 4

# conversion of the fits files. Should only be run once.
do_convert_fits=0
# perform photometry on converted fits files. Should only be run once.
do_photometry=0
# perform matching of photometry files. Should only be run once.
do_match=0
# search for ideal aperture for these matched photometry files. Results are saved when run, so can be set to 0 as long
# as photometry/matching doesn't change
do_aperture_search=1
# write X/Y and RA/DEC positions of each star to disk
do_pos=0
# process compstars
do_compstars=1
# write lightcurves to disk
do_lightcurve=1
# perform machine learning to find new variables
do_ml=0
# write lightcurve plots
do_lightcurve_plot=0
# write phase diagrams
do_phase_diagram=0
# write various other charts
do_field_charts=0
# perform aavso reporting
do_reporting=1

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
testdir = basedir + 'test/'
conf_phot=basedir+'muniphot.conf'
conf_match=basedir+'match.conf'

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

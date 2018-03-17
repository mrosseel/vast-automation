#### DO NOT TOUCH ####
import multiprocessing
codedir = './'

aperture = 2
all = range(1,10000)
thou = range(1,1000)
selection = [1,2,3,138]
selection2 = (1,2,3,73,138,143,264,2675,1045,847,1193)
custom_muniwin = range(5776,10000)
#custom_charts = range(1,5776)
custom_charts = range(1,2164)
all_star_list = range(1,138)
#### DO NOT TOUCH ####

#### SETTINGS for DO_MUNIWIN.PY ####

#basedir = codedir + 'inputfiles/WWCrA_allflat/'
basedir =  codedir + 'inputfiles/WWCrA2015/'
star_list = thou
nr_threads = multiprocessing.cpu_count() * 2

do_convert_fits=False
do_photometry=False
do_match=False
do_munifind=False
do_lightcurve=False
do_lightcurve_resume=False
do_pos=False
do_pos_resume=False
do_calibrate=False
do_upsilon=False
do_naming=False
do_charting=False
do_phase_diagram=True

### CALIBRATION ###
reference_frame = 'WWCrA#30V_000184527_FLAT.fit'
reference_header = basedir + 'new-image.fits'
ra_deg = 271.4032917
dec_deg = -43.8326111
### CALIBRATION

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
#### DO NOT TOUCH ####

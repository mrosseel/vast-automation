import os
import shutil

#basedir = os.getcwd()+'/inputfiles/WWCrA_allflat/'
basedir = os.getcwd()+'/inputfiles/WWCrA2015/'
fitsdir = basedir + "fits/"
convfitsdir = basedir + "converted_fits/"
photometrydir = basedir + "photometry/"
matchedphotometrydir = basedir + "matched_photometry/"
resultdir = basedir + "results/"
lightcurvedir = resultdir + "lightcurves/"
posdir = resultdir + "positions/"

### CALIBRATION ###
reference_frame = basedir + 'WWCrA#30V_000184527_FLAT.fit'
xpos = 717
ypos = 654
reference_header = basedir + 'wcs.fits'
reference_star = 73
arcsecond_x = 47*60
arcsecond_y = 47*60
### CALIBRATION

aperture = 2
all = range(1,10000)
thou = range(1,1000)
selection = [1,2,3,137]
custom_muniwin = range(5776,10000)
#custom_charts = range(1,5776)
custom_charts = range(1,2164)
all_star_list = range(1,138)

def trash_and_recreate_dir(dir):
    shutil.rmtree(dir, ignore_errors=True)
    os.makedirs(dir, exist_ok=True)

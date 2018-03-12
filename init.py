codedir = './'

#basedir = codedir + 'inputfiles/WWCrA_allflat/'
basedir =  codedir + 'inputfiles/WWCrA2015/'

match_with_photometry_file = 'phot000001.pht'
do_match=False
do_munifind=False
do_lightcurve=True
do_pos=True
do_pos_resume=False
do_calibrate=True
do_charts=False

### CALIBRATION ###
reference_frame = 'WWCrA#30V_000184527_FLAT.fit'
reference_header = basedir + 'new-image.fits'
xpos = 962.856994023
ypos = 711.327534113
ra_deg = 271.4032917
dec_deg = -43.8326111
### CALIBRATION

aperture = 2
all = range(1,10000)
thou = range(1,1000)
selection = [1,2,3,137]
selection2 = (1,2,3,73,137,143,264,2675,1045,847,1193)
custom_muniwin = range(5776,10000)
#custom_charts = range(1,5776)
custom_charts = range(1,2164)
all_star_list = range(1,138)
star_list = thou

# standard directories
fitsdir = basedir + "fits/"
convfitsdir = basedir + "converted_fits/"
photometrydir = basedir + "photometry/"
matchedphotometrydir = basedir + "matched_photometry/"
resultdir = basedir + "results/"
lightcurvedir = resultdir + "lightcurves/"
posdir = resultdir + "positions/"
worldposdir = resultdir + "world_positions/"

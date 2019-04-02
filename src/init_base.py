#### DO NOT TOUCH ####
# standard directories
datadir = 'current/'
reference_file = 'new-image.fits'
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

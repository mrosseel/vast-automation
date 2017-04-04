import os

basedir = os.getcwd()+'/inputfiles/WWCrA/'
lightcurve_dir = basedir + "outstars/"
#reference_frame = basedir+'WWCrA#30V_000409679_FLAT.fit'
reference_frame = 'wcs.fits' # this is a fits file calculated by astrometry.net
aperture = 2
one_tenthousand = range(1,10000)
two = [1,2]
all_star_list = two

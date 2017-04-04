import os

basedir = os.getcwd()+'/inputfiles/WWCrA2015/'
lightcurve_dir = basedir + "outstars/"
#reference_frame = basedir+'WWCrA#30V_000409679_FLAT.fit'
reference_frame = 'wcs.fits' # this is a fits file calculated by astrometry.net
aperture = 2
all = range(1,10000)
thou = range(1,1000)
two = [1,2]
all_star_list = thou

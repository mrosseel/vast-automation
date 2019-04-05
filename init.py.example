from multiprocessing import cpu_count

all = range(1, 10001)
thousand = range(1, 1000)
hundred = range(1, 100)

#### SETTINGS for DO_MUNIWIN.PY ####


aperture_find_percentage=0.05
free_memory_GB=10
star_list = hundred
nr_threads = cpu_count()

# Comparison stars to use:
comparison_stars = ['UCAC4 231-154752', 'UCAC4 232-147677']
# use None if you want automatic aperture
aperture = 4

# conversion of the fits files. Should only be run once.
do_convert_fits=1
# perform photometry on converted fits files. Should only be run once.
do_photometry=1
# perform matching of photometry files. Should only be run once.
do_match=1
# search for ideal aperture for these matched photometry files. Results are saved when run, so can be set to 0 as long
# as photometry/matching doesn't change
do_aperture_search=1
# write X/Y and RA/DEC positions of each star to disk
do_pos=1
# process compstars
do_compstars=1
# write lightcurves to disk
do_lightcurve=1
# perform machine learning to find new variables
do_ml=1
# write lightcurve plots
do_lightcurve_plot=1
# write phase diagrams
do_phase_diagram=1
# write various other charts
do_field_charts=1
# perform aavso reporting
do_reporting=1

### CALIBRATION ###

ra_deg = 271.4032917
dec_deg = -43.8326111
sitelat = '-22 57 10'
sitelong = '-68 10 49'
sitealt= 2500
### CALIBRATION



import os

basedir = os.getcwd()+'/inputfiles/WWCrA_allflat/'
lightcurve_dir = basedir + "outstars/"
#reference_frame = basedir+'WWCrA#30V_000409679_FLAT.fit'
reference_frame = 'wcs.fits' # this is a fits file calculated by astrometry.net
aperture = 2
all = range(1,10000)
thou = range(1,1000)
selection = [1,2,3,137]
custom_muniwin = range(5776,10000)
#custom_charts = range(1,5776)
custom_charts = range(1,2164)
all_star_list = [73, 2411, 2804, 4227, 5155, 797, 2337, 2587, 4074, 5506, 1103, 3722, 2976, 2539, 142, 452, 1520, 5256, 316, 4944, 4479, 721, 4346, 3516, 2829, 5391, 5545, 4247, 268, 5572, 1473, 2757, 3879, 5720, 1930, 3523, 4835, 1489, 4898, 5590, 1759, 1780, 1804, 3394, 2120, 3578, 3913, 886, 2260, 3007, 4676, 4920, 5037, 5671, 2519, 3488, 5580, 3161, 5007, 5187, 344, 354, 743, 849, 2878, 2985, 3782, 4244, 4852, 4959, 5022, 5133, 5224, 1067, 1495, 2126, 2626, 2665, 3265, 3718, 3934, 5114, 5277, 1610, 3340, 3610, 3957, 4156, 4725, 5523, 5, 773, 3290, 3681, 3954, 4912, 5271, 5526, 5687]

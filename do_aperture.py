import init
from do_muniwin import getBestComparisonStars
from do_muniwin import write_munifind
from reading import trash_and_recreate_dir
import numpy as np
import pandas as pd

# match_file = 'match000???.pht'
def find_optimal_aperture(match_file, aperture_range = init.aperture_range):
    trash_and_recreate_dir(init.aperturedir)
    print('Finding optimal aperture for match file {}, range {}'.format(match_file, aperture_range))
    results = []
    for aperture in aperture_range:
        print('Calculating aperture:', aperture)
        current_file = write_munifind(aperture, match_file, quiet=True)
        try:
            print("Getting best comparison stars from ", current_file)
            result = getBestComparisonStars(100, current_file)
            results.append([aperture, result['STDEV'].mean()])
        except:
            results.append([aperture, 100])
    diff_result_array = []
    for mean in results:
        diff_result_array = np.append(diff_result_array, mean)
    optimal_aperture = aperture_range[pd.Series(diff_result_array).idxmin()]
    print('Optimal aperture:', optimal_aperture)
    np.savetxt(init.aperturedir + "aperture.csv", diff_result_array, delimiter=",")
    np.savetxt(init.aperturedir + "aperture_best.txt", [optimal_aperture])
    return optimal_aperture

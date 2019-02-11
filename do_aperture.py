import init
from do_muniwin import getBestComparisonStars
from do_muniwin import write_munifind
from reading import trash_and_recreate_dir
import numpy as np
import pandas as pd
import tqdm
import multiprocessing as mp
from functools import partial

# match_file = 'match000???.pht'
def find_optimal_aperture(match_file, aperture_range = init.aperture_range):
    trash_and_recreate_dir(init.aperturedir)
    print('Finding optimal aperture for match file {}, range {}'.format(match_file, aperture_range))
    # rewrite using https://docs.python.org/3/library/concurrent.futures.html
    pool = mp.Pool(1, maxtasksperchild=100)
    results = []
    func = partial(calculate_aperture, match_file=match_file)
    for resulttuple in tqdm.tqdm(pool.imap_unordered(func, aperture_range, 1), total=len(aperture_range)):
        results.append(resulttuple)
    results = sorted(results, key=lambda x: x[0])
    print(results)
    diff_result_array = [x[1] for x in results]
    print(diff_result_array)
    print(f"diff_result_array: {len(diff_result_array)}, aperture_range: {len(aperture_range)}")
    optimal_aperture = aperture_range[pd.Series(diff_result_array).idxmin()]
    print('Optimal aperture:', optimal_aperture)
    np.savetxt(init.aperturedir + "aperture.csv", diff_result_array, delimiter=",")
    np.savetxt(init.aperturedir + "aperture_best.txt", [optimal_aperture])
    return optimal_aperture

# match_file = 'match000???.pht'
def calculate_aperture(aperture, match_file):
        print('Calculating aperture:', aperture)
        current_file = write_munifind(aperture, match_file=True, quiet=True, percentage=init.munifind_percentage)
        try:
            result = getBestComparisonStars(100, current_file)
            return [aperture, result['STDEV'].mean()]
        except:
            return [aperture, 100]

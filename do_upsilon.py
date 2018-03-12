import upsilon_helper
import init
from multiprocessing import Pool, Queue
import tqdm

from functools import partial

def start_upsilon(star_list, star_limit):
    #print("Starting upsilon with nr stars", len(star_list))
    pool = Pool(init.nr_threads)
    func = partial(upsilon_helper.predict_star,limit=star_limit)
    result_list = []
    #print("Predicting variability for ",len(star_list),"stars")
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list, chunksize=10), total=len(star_list)):
        result_list.append(_)
        pass
    return result_list

def run(star_list):
    upsilon_output =init.basedir+'upsilon_output.txt'
    unsorted_result_list = start_upsilon(star_list, 5000)
    # save before sorting
    upsilon_helper.save_results(unsorted_result_list, upsilon_output)
    # sort on probability, descending - primary key
    sorted_result_list = sorted(unsorted_result_list, key=lambda result: result[2], reverse=True)
    # sort on danger flag, ascending - secondary key
    sorted_result_list = sorted(sorted_result_list, key=lambda result: result[3])
    # save #2 after sorting
    upsilon_helper.save_results(sorted_result_list, upsilon_output)

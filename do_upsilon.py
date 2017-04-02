import upsilon_helper
import init
import upsilon
from multiprocessing import Pool, Queue
import tqdm

from functools import partial

def start_upsilon(star_list):
    print("Starting upsilon with nr stars", len(star_list))
    pool = Pool(8)
    func = partial(upsilon_helper.predict_star)
    result_list = []
    print("Predicting variability for ",len(star_list),"stars")
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list, chunksize=10), total=len(star_list)):
        result_list.append(_)
        pass
    return result_list

my_result_list = start_upsilon(range(1,100))
upsilon_output =init.basedir+'upsilon_output.txt'
upsilon_helper.save_results(my_result_list, upsilon_output)

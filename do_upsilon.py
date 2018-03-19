import init
import upsilon
import reading
from multiprocessing import Pool, Queue
import tqdm
import pandas as pd
from functools import partial

def start_upsilon(star_list, star_limit):
    #print("Starting upsilon with nr stars", len(star_list))
    pool = Pool(init.nr_threads)
    func = partial(predict_star,limit=star_limit)
    result_list = []
    #print("Predicting variability for ",len(star_list),"stars")
    for _ in tqdm.tqdm(pool.imap_unordered(func, star_list, chunksize=10), total=len(star_list)):
        result_list.append(_)
        pass
    return result_list

def test_upsilon():
    upsilon.test_predict()

def save_results(result_list, output_file):
    final_list = []
    column_names = ['star', 'label', 'probability', 'flag']
    columns_done = False
    for entry in result_list:
        features_dict = entry[4]
        entry = entry[:-1]
        for key in features_dict:
            entry.append(features_dict[key])
            if not columns_done: column_names.append(key)
        final_list.append(entry)
        columns_done = True
    df=pd.DataFrame(final_list,columns=column_names)
    print(df.head())
    df.to_csv(output_file, index=False)

# returns [ star, label, probability, flag ]
def predict_star(star, limit=-1):
    #print("star:",star)
    try:
        df = reading.read_lightcurve(star)
        if(limit > 0):
            df = df[:limit]
            print("Restricting star", star, " to limit:", limit)
        mag = df['V-C']
        date = df['JD']
        err = df['s1']
        e_features = upsilon.ExtractFeatures(date, mag, err)
        e_features.run()
        features = e_features.get_features()

        # Classify the light curve
        label, probability, flag = upsilon.predict(rf_model, features)
        return [star, label, probability, flag, features]
    except:
        print(star, 'error')
        return [star, 'NA', -1, 1]

def run(star_list):
    upsilon_output =init.basedir+'upsilon_output.txt'
    unsorted_result_list = start_upsilon(star_list, 5000)
    # save before sorting
    save_results(unsorted_result_list, upsilon_output)
    # sort on probability, descending - primary key
    sorted_result_list = sorted(unsorted_result_list, key=lambda result: result[2], reverse=True)
    # sort on danger flag, ascending - secondary key
    sorted_result_list = sorted(sorted_result_list, key=lambda result: result[3])
    # save #2 after sorting
    save_results(sorted_result_list, upsilon_output)

rf_model = upsilon.load_rf_model()

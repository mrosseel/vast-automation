#cell 16
import init
import upsilon
import pandas as pd
import reading

def test_upsilon():
    upsilon.test_predict()

#no longer used
def predict_star_list(all_star_list, output_file):
    rf_model = upsilon.load_rf_model()
    full_list = []
    for star in all_star_list:
        try:
            df = reading.read_lightcurve(star)
            mag = df['V-C']
            date = df.index.values
            e_features = upsilon.ExtractFeatures(date, mag)
            e_features.run()
            features = e_features.get_features()

            # Classify the light curve
            label, probability, flag = upsilon.predict(rf_model, features)
            #print("star, label, probability, flag")
            #print(star, label, probability, flag)
            full_list.append([star, label, probability, flag])
            #print("result list", my_result_list)
        except:
            print(star, 'error')
            full_list.append([star, 'NA', 'NA', 'NA'])

def save_results(result_list, output_file):
    df=pd.DataFrame(result_list,columns=['star', 'label', 'probability', 'flag'])
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
        #date = df.index.values
        date = df['JD']
        err = df['s1']
        e_features = upsilon.ExtractFeatures(date, mag, err)
        e_features.run()
        features = e_features.get_features()
        #print(features)

        # Classify the light curve
        label, probability, flag = upsilon.predict(rf_model, features)
        #print("star, label, probability, flag")
        #print(star, label, probability, flag)
        return [star, label, probability, flag, features]
    except:
        print(star, 'error')
        return [star, 'NA', -1, 1]

rf_model = upsilon.load_rf_model()

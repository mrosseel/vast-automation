#cell 16
import init
import upsilon
import pandas as pd

def test_upsilon():
    upsilon.test_predict()

def read_lightcurve(star):
    #print("Reading lightcurve", star)
    df = pd.read_csv(init.lightcurve_dir + str(star) + '.txt', skiprows=[1], sep=' ')
    df = df[df['V-C'] < 99]
    return df

def read_pos(star):
    df = pd.read_csv(init.lightcurve_dir + 'pos_' + str(star) + '.txt', skiprows=[1], sep=' ')
    df = df[df['X'] > 0]
    df = df[df['MAG'] < 99]
    return (df['X'].iloc[1], df['Y'].iloc[1])

def predict_star_list(all_star_list, output_file):
    rf_model = upsilon.load_rf_model()
    full_list = []
    for star in all_star_list:
        try:
            df = read_lightcurve(star)
            mag = df['V-C']
            date = df.index.values
            e_features = upsilon.ExtractFeatures(date, mag)
            e_features.run()
            features = e_features.get_features()

            # Classify the light curve
            label, probability, flag = upsilon.predict(rf_model, features)
            print("star, label, probability, flag")
            print(star, label, probability, flag)
            full_list.append([star, label, probability, flag])
            print("result list", my_result_list)
        except:
            print(star, 'error')
            full_list.append([star, 'NA', 'NA', 'NA'])

def save_results(result_list, output_file):
    df=pd.DataFrame(result_list,columns=['star', 'label', 'probability', 'flag'])
    print(df.head())
    df.to_csv(output_file, index=False)

def predict_star(star):
    print("star:",star)
    try:
        df = read_lightcurve(star)
        mag = df['V-C']
        date = df.index.values
        e_features = upsilon.ExtractFeatures(date, mag)
        e_features.run()
        features = e_features.get_features()

        # Classify the light curve
        label, probability, flag = upsilon.predict(rf_model, features)
        #print("star, label, probability, flag")
        #print(star, label, probability, flag)
        return [star, label, probability, flag]
    except:
        print(star, 'error')
        return [star, 'NA', 'NA', 'NA']

rf_model = upsilon.load_rf_model()

import pickle
import pandas as pd
import argparse
import tqdm


def convert(path):
    print(f"Reading {path}...")
    # using pandas with a column specification
    col_specification =[(0, 7), (8, 38), (39, 40), (41,50), (51,60), (61,91), (93,94), (94,100),(100, 101), (101,107), (108,109), (109,110), (110,116), (116,117), (117,123), (124,136), (136,137), (138,139), (139,155), (155,158)]
    data = pd.read_fwf(path, colspecs=col_specification, names=('OID', 'Name','V','RAdeg', 'DEdeg', 'Type', 'l_max', 'max', 'u_max', 'n_max', 'f_min', 'l_min', 'min', 'u_min', 'n_min', 'Epoch', 'u_Epoch', 'l_Period', 'Period', 'u_Period'), converters={'l_max': str})

    ra_deg_np = data['RAdeg'].values
    dec_deg_np = data['DEdeg'].values
    metadata = []
    for index, row in tqdm.tqdm(data.iterrows(), total=len(data), desc="Converting"):
        metadata.append({'id': index, 'OID': row['OID'], 'Name': row['Name'], 'Type': row['Type'], 'l_Period': row['l_Period'], 'Period': row['Period'], 'u_Period': row['u_Period']})
    vsx_dict = { 'ra_deg_np': ra_deg_np, 'dec_deg_np': dec_deg_np, 'metadata': metadata}

    outfile = './vsx_catalog.bin'
    print(f"Writing {outfile}...")
    # From Python 3.6 onwards, the standard dict type maintains insertion order by default. => no OrderedDict necessary
    with open(outfile, 'wb') as fp:
        pickle.dump(vsx_dict, fp)
    print(f"Done.")


# returns:   { 'ra_deg_np': ra_deg_np, 'dec_deg_np': dec_deg_np, 'metadata': ({'id': index, 'OID': row['OID'],
# 'Name': row['Name'], 'Type': row['Type'], 'l_Period': row['l_Period'], 'Period': row['Period'],
# 'u_Period': row['u_Period']})}
def read(path):
    with open(path, 'rb') as fp:
        read_dict = pickle.load(fp, encoding='latin1')
    return read_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converting vsx to pickle file vsx.bin')
    parser.add_argument('vsx_path')
    args = parser.parse_args()
    convert(args.vsx_path)

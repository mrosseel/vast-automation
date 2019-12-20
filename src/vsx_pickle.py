import pickle
import pandas as pd
import argparse
import tqdm


# http://cdsarc.u-strasbg.fr/ftp/B/vsx/ReadMe

# Byte-by-byte Description of file: vsx.dat
# --------------------------------------------------------------------------------
# Bytes Format Units   Label    Explanations
# --------------------------------------------------------------------------------
# 1-  7  I7    ---     OID      Internal identifier, can be used to
# link out to the VSX database (1)
# 9- 38  A30   ---     Name     Variable star identifier
# 40  I1    ---     V        [0,3] Variability flag (2)
# 42- 50  F9.5  deg     RAdeg    Right ascension (J2000)
# 52- 60  F9.5  deg     DEdeg    Declination (J2000)
# 62- 91  A30   ---     Type     Variability type, as in GCVS catalog
# 94  A1    ---   l_max      Limit flag on max
# 95-100  F6.3  mag     max      ? Magnitude at maximum, or amplitude
# 101  A1    ---   u_max      Uncertainty flag on max
# 102-107  A6    ---   n_max      Passband on max magnitude (4)
# 109  A1    ---   f_min      [(] '(' indicates an amplitude
#                               110  A1    ---   l_min      Limit flag on min
#                               111-116  F6.3  mag     min      ? Magnitude at minimum, or amplitude
#                               117  A1    ---   u_min      Uncertainty flag on min
#                               118-123  A6    ---   n_min      Passband on min magnitude (4)
# 125-136  F12.4 d       Epoch    ? Epoch of maximum or minimum (HJD)
# 137  A1    ---   u_Epoch    [:)] Uncertainty flag (:) on epoch
# 139  A1    ---   l_Period   [<>(] Limit flag on period (3)
# 140-155 F16.10 d       Period   ? Period of the variable in days
# 156-158  A3    ---   u_Period   [:)*/N2 ] Uncertainty flag on Period (3)
# --------------------------------------------------------------------------------

def convert(args):
    path = args.vsx_path
    outfile = './vsx_catalog.bin'
    print(f"Reading {path}...")
    # using pandas with a column specification
    col_specification = [(0, 7), (8, 38), (39, 40), (41, 50), (51, 60), (61, 91), (93, 94), (94, 100), (100, 101),
                         (101, 107), (108, 109), (109, 110), (110, 116), (116, 117), (117, 123), (124, 136), (136, 137),
                         (138, 139), (139, 155), (155, 158)]
    data = pd.read_fwf(path, colspecs=col_specification, skiprows=28, names=(
    'OID', 'Name', 'V', 'RAdeg', 'DEdeg', 'Type', 'l_max', 'max', 'u_max', 'n_max', 'f_min', 'l_min', 'min', 'u_min',
    'n_min', 'Epoch', 'u_Epoch', 'l_Period', 'Period', 'u_Period'), converters={'l_max': str})

    if args.test:
        data = data.loc[0:10]
        outfile = './vsx_mini.bin'
    ra_deg_np = data['RAdeg'].values
    dec_deg_np = data['DEdeg'].values
    extradata = data.to_dict(orient='records')
    vsx_dict = {'ra_deg_np': ra_deg_np, 'dec_deg_np': dec_deg_np, 'extradata': extradata}

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
    parser.add_argument('-t', '--test', help="Create a mini test vsx file", required=False, action='store_true')
    parser.add_argument('vsx_path')
    args = parser.parse_args()
    convert(args)

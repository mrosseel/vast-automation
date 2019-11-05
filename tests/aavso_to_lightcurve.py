import pandas as pd
import utils

def run():
    df = pd.read_csv('./data/outlier.aavso', sep=',', skiprows=7, header=None, index_col=False,
                     names=['NAME', 'DATE', 'MAG', 'MERR', 'FILT', 'TRANS', 'MTYPE', 'CNAME', 'CMAG', 'KNAME', 'KMAG', 'AMASS',
                            'GROUP', 'CHART,NOTES'])
    print(df.head())
    print(df['MAG'].head())
    df2 = utils.reject_outliers_iqr(df, 'MAG', 5)
    print(len(df), len(df2))


if __name__ == '__main__':
    run()

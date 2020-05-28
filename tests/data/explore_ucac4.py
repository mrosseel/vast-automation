import ucac4
from pathlib import Path
import logging
import pandas as pd
from tqdm.auto import trange

def process_file(uc, zone_nr, bucket, process_func):
    filecontent = uc.get_zone_filecontent(zone_nr)
    stars = []
    logging.debug(f"zone: {zone_nr}, bucket: {bucket}")
    run_nrs = uc.index_bin_to_run_nrs(zone_nr, bucket)
    logging.debug(f"run nrs: {run_nrs}")
    if len(run_nrs) == 0:
        logging.debug(f"empty run {zone_nr}, {bucket}, {run_nrs}")
        return
    for run_nr in run_nrs:
        # logging.debug(f"unpack: {uc.zone_star_length * run_nr}:{uc.zone_star_length * (run_nr+1) +1}")
        result = uc.unpack_zone_fileformat(filecontent[uc.zone_star_length * run_nr:uc.zone_star_length * (run_nr+1)])
        star = uc.make_startuple(result)
        stars.append(star)
    return process_func(uc, zone_nr, bucket, run_nrs, stars)


def print_pretty_stars(uc, zone_nr, bucket, run_nrs, stars):
    for star in stars:



def get_ra_decs(uc, zone_nr, bucket, run_nrs, stars):
    ################ DEBUG
    radecs = [uc.get_real_ra_dec(x.ra, x.spd) for x in stars]
    rasraw = [x.ra for x in stars]
    decsraw = [x.spd for x in stars]
    ras = [x[0] for x in radecs]
    decs = [x[1] for x in radecs]
    logging.debug(f"Zone {zone_nr}: between ra<{min(ras)}, {max(ras)}> and  dec<f{min(decs)}, {max(decs)}>, "
                  f"raw: ra<{min(rasraw)}, {max(rasraw)}> and  dec<{min(decsraw)}, {max(decsraw)}>")
    logging.debug(f"run nrs: {run_nrs}")
    return [zone_nr, bucket, min(ras), max(ras), min(decs), max(decs), min(rasraw), max(rasraw), min(decsraw),
            max(decsraw)]


def do_single():
    uc = ucac4.UCAC4(ucac_path=Path('../support/ucac4/UCAC4/'))
    create_df([process_file(uc, 7, 1100)])


def do_zone():
    uc = ucac4.UCAC4(ucac_path=Path('../support/ucac4/UCAC4/'))
    result = []
    for bucket in range(0, 1439):
        result.append(process_file(uc, 231, bucket))
    create_df(result)

def do_bucket():
    uc = ucac4.UCAC4(ucac_path=Path('../support/ucac4/UCAC4/'))
    result = []
    for zone in range(1, 901):
        result.append(process_file(uc, zone, 1))
    create_df(result)


def do_all():
    uc = ucac4.UCAC4(ucac_path=Path('../support/ucac4/UCAC4/'))
    result = []
    for count in trange(1, 900, desc='zones'):
        for run_nr in trange(1, 1440, desc='buckets', leave=False):
            result.append(process_file(uc, count, run_nr, get_ra_decs()))
    create_df(result)


def create_df(resultarr):
    resultarr = [x for x in resultarr if x is not None]
    # print(resultarr)
    df = pd.DataFrame(resultarr, columns=['zone_nr', 'bucket', 'min_ras', 'max_ras', 'min_decs', 'max_decs',
                                          'min_rasraw', 'max_rasraw', 'min_decsraw', 'max_decsraw'])
    df = df.dropna()
    df.to_csv('zonebucket.csv')
    print("done")
    # for entry in range(1, 1440):
    #     print(f" {df[]}")
    # print(df.describe(), df.info())
    # print(df['min_ras'].min(), df['max_ras'].max())


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.WARN)
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s")
    # uc = ucac4.UCAC4(ucac_path=Path('../support/ucac4/UCAC4/'))
    do_all()

# col byte item   fmt unit       explanation                            notes
# ---------------------------------------------------------------------------
# 1  1- 3 ra     I*4 mas        right ascension at  epoch J2000.0 (ICRS) (1)
# 2  5- 8 spd    I*4 mas        south pole distance epoch J2000.0 (ICRS) (1)
# 3  9-10 magm   I*2 millimag   UCAC fit model magnitude                 (2)
# 4 11-12 maga   I*2 millimag   UCAC aperture  magnitude                 (2)
# 5 13    sigmag I*1 1/100 mag  error of UCAC magnitude                  (3)
# 6 14    objt   I*1            object type                              (4)
# 7 15    cdf    I*1            combined double star flag                (5)
# 15 bytes
# 8 16    sigra  I*1 mas        s.e. at central epoch in RA (*cos Dec)   (6)
# 9 17    sigdc  I*1 mas        s.e. at central epoch in Dec             (6)
# 10 18    na1    I*1            total # of CCD images of this star
# 11 19    nu1    I*1            # of CCD images used for this star       (7)
# 12 20    cu1    I*1            # catalogs (epochs) used for proper motions
# 5 bytes
# 13 21-22 cepra  I*2 0.01 yr    central epoch for mean RA, minus 1900
# 14 23-24 cepdc  I*2 0.01 yr    central epoch for mean Dec,minus 1900
# 15 25-26 pmrac  I*2 0.1 mas/yr proper motion in RA*cos(Dec)             (8)
# 16 27-28 pmdc   I*2 0.1 mas/yr proper motion in Dec
# 17 29    sigpmr I*1 0.1 mas/yr s.e. of pmRA * cos Dec                   (9)
# 18 30    sigpmd I*1 0.1 mas/yr s.e. of pmDec                            (9)
# 10 bytes
# 19 31-34 pts_key I*4           2MASS unique star identifier            (10)
# 20 35-36 j_m    I*2 millimag   2MASS J  magnitude
# 21 37-38 h_m    I*2 millimag   2MASS H  magnitude
# 22 39-40 k_m    I*2 millimag   2MASS K_s magnitude
# 23 41    icqflg I*1            2MASS cc_flg*10 + ph_qual flag for J    (11)
#     24 42     (2)   I*1            2MASS cc_flg*10 + ph_qual flag for H    (11)
#     25 43     (3)   I*1            2MASS cc_flg*10 + ph_qual flag for K_s  (11)
#     26 44    e2mpho I*1 1/100 mag  error 2MASS J   magnitude               (12)
# 27 45     (2)   I*1 1/100 mag  error 2MASS H   magnitude               (12)
# 28 46     (3)   I*1 1/100 mag  error 2MASS K_s magnitude               (12)
# 16 bytes
# 29 47-48 apasm  I*2 millimag   B magnitude from APASS                  (13)
# 30 49-50  (2)   I*2 millimag   V magnitude from APASS                  (13)
# 31 51-52  (3)   I*2 millimag   g magnitude from APASS                  (13)
# 32 53-54  (4)   I*2 millimag   r magnitude from APASS                  (13)
# 33 55-56  (5)   I*2 millimag   i magnitude from APASS                  (13)
# 34 57    apase  I*1 1/100 mag  error of B magnitude from APASS         (14)
# 35 58     (2)   I*1 1/100 mag  error of V magnitude from APASS         (14)
# 36 59     (3)   I*1 1/100 mag  error of g magnitude from APASS         (14)
# 37 60     (4)   I*1 1/100 mag  error of r magnitude from APASS         (14)
# 38 61     (5)   I*1 1/100 mag  error of i magnitude from APASS         (14)
# 39 62    gcflg  I*1            Yale SPM g-flag*10  c-flag              (15)
# 16 bytes
# 40 63-66 icf(1) I*4            FK6-Hipparcos-Tycho source flag         (16)
# 41       icf(2) ..             AC2000       catalog match flag         (17)
# 42       icf(3) ..             AGK2 Bonn    catalog match flag         (17)
# 43       icf(4) ..             AKG2 Hamburg catalog match flag         (17)
# 44       icf(5) ..             Zone Astrog. catalog match flag         (17)
# 45       icf(6) ..             Black Birch  catalog match flag         (17)
# 46       icf(7) ..             Lick Astrog. catalog match flag         (17)
# 47       icf(8) ..             NPM  Lick    catalog match flag         (17)
# 48       icf(9) ..             SPM  YSJ1    catalog match flag         (17)
# 4 bytes
# 49 67    leda   I*1            LEDA galaxy match flag                  (18)
# 50 68    x2m    I*1            2MASS extend.source flag                (19)
# 51 69-72 rnm    I*4            unique star identification number       (20)
# 52 73-74 zn2    I*2            zone number of UCAC2 (0 = no match)     (21)
# 53 75-78 rn2    I*4            running record number along UCAC2 zone  (21)
# 12 bytes
# ---------------------------------------------------------------------------
# 78 = total number of bytes per star record
#

# C code converted to correct prefixes:

# i ra
# i spd;         /* RA/dec at J2000.0,  ICRS,  in milliarcsec */
# H mag1
# H mag2;     /* UCAC fit model & aperture mags, .001 mag */
# B mag_sigma;
# B obj_type
# B double_star_flag;
#
# b ra_sigma
# b dec_sigma;    /* sigmas in RA and dec at central epoch */
# B n_ucac_total;      /* Number of UCAC observations of this star */
# B n_ucac_used;      /* # UCAC observations _used_ for this star */
# B n_cats_used;      /* # catalogs (epochs) used for prop motion */
#
# H epoch_ra;        /* Central epoch for mean RA, minus 1900, .01y */
# H epoch_dec;       /* Central epoch for mean DE, minus 1900, .01y */
# h pm_ra;            /* prop motion, .1 mas/yr = .01 arcsec/cy */
# h pm_dec;           /* prop motion, .1 mas/yr = .01 arcsec/cy */
# b pm_ra_sigma;       /* sigma in same units */
# b pm_dec_sigma;
#
# I twomass_id;        /* 2MASS pts_key star identifier */
# H mag_j
# H mag_h
# H mag_k;  /* 2MASS J, H, K_s mags,  in millimags */
# BBB icq_flag[3];
# BBB e2mpho[3];          /* 2MASS error photometry (in centimags) */
# HHHHH apass_mag[5];      /* in millimags */
# bbbbb apass_mag_sigma[5];  /* in centimags */
# B yale_gc_flags;      /* Yale SPM g-flag * 10 + c-flag */
# I catalog_flags;
# B leda_flag;          /* LEDA galaxy match flag */
# B twomass_ext_flag;   /* 2MASS extended source flag */
# I id_number;
# H ucac2_zone;
# I ucac2_number;

# NOTE !!!!!!!!!!!  The ra/dec and proper motion sigmas are now offset by 128; i.e.,  a value of -128 would indicate a zero sigma.
import logging
import struct
import math
import numpy as np
from collections import namedtuple
from typing import List, Tuple
from star_description import StarDescription
from astropy.coordinates import SkyCoord
from pathlib import Path
from LRUCache import LRUCache
# remove this, would be cleaner in ucac4_utils or something
import do_calibration
import tqdm

StarTuple = namedtuple('Star', 'ra spd mag1 mag2 mag_sigma obj_type double_star_flag ra_sigma'
                               ' dec_sigma n_ucac_total n_ucac_used n_cats_used epoch_ra epoch_dec pm_ra pm_dec'
                               ' pm_ra_sigma pm_dec_sigma twomass_id mag_j mag_h mag_k icq_flag1 icq_flag2 icq_flag3'
                               ' e2mpho1 e2mpho2 e2mpho3 apass_mag_B apass_mag_V apass_mag_g apass_mag_r apass_mag_i'
                               ' apass_mag_sigma_B apass_mag_sigma_V apass_mag_sigma_g apass_mag_sigma_r'
                               ' apass_mag_sigma_i yale_gc_flags catalog_flags leda_flag twomass_ext_flag'
                               ' id_number ucac2_zone ucac2_number')


def get_line_nr(n0, nn, line):
    return f"{n0[line * 900]}\t{nn[line * 900]}"


class UCAC4:
    def __init__(self, ucac_path: Path = Path('./support/ucac4/UCAC4/')):
        # star id given by munipack
        self.ucac_path = ucac_path
        self.zones_cache = LRUCache(capacity=5)
        self.bucket_cache = LRUCache(capacity=100000)
        self.zone_bucket_cache = LRUCache(capacity=10000)
        self.index_cache = None
        with open(str(Path(ucac_path, 'u4i', 'u4index.unf')), mode='rb') as file:  # b is important -> binary
            self.index_cache = file.read()
        self.zone_starformat = "=iiHHBBBbbBBBHHhhbbIHHHBBBBBBHHHHHbbbbbBIBBIHI"
        self.zone_star_length = struct.calcsize(self.zone_starformat)
        self.unpack_zone_fileformat = struct.Struct(self.zone_starformat).unpack
        self.result_n0, self.result_nn = UCAC4._get_n0_and_nn(self.index_cache)


    def get_ucac4_details(self, ucac_id) -> List[Tuple[Tuple, str, int]]:
        full_id = ucac_id[6:]
        zone = full_id[:3]
        run_nr = int(full_id[4:].lstrip("0"))
        logging.debug(f"UCAC4 id {full_id}, {zone}, {run_nr}")
        return self.get_ucac4_details_raw(zone, [run_nr])


    def get_zone_from_file(self, zone: str):
        result = self.zones_cache.get(zone)
        if result != -1:
            return result
        with open(str(Path(self.ucac_path, f"u4b/z{zone:03}")), mode='rb') as file:  # b is important -> binary
            result = file.read()
            self.zones_cache.set(zone, result)
            return result


    def get_ucac4_details_raw(self, zone: str, run_nr: List[int]) -> List[Tuple[StarTuple, str, int]]:
        stars = []
        filecontent = self.get_zone_from_file(zone)
        for number in run_nr:
            key = (zone, number)
            bucketcache = self.bucket_cache.get(key)
            if bucketcache != -1:
                star = bucketcache
            else:
                result = self.unpack_zone_fileformat(
                    filecontent[self.zone_star_length * (number - 1):self.zone_star_length * number])
                star = UCAC4.make_star(result)
                self.bucket_cache.set(key, star)
            stars.append((star, zone, number))
        # logging.debug(f"read ucac4 star: {star}")
        return stars

    @staticmethod
    def _get_n0_and_nn(index_cache):
        index_format = "1296000I"
        index_length = struct.calcsize(index_format)
        logging.debug(f"index length  is {index_length}")
        unpack_starformat = struct.Struct(index_format).unpack
        result_n0 = unpack_starformat(index_cache[0:index_length])
        result_nn = unpack_starformat(index_cache[index_length:])
        return result_n0, result_nn


    def get_ucac4_sd(self, ra: float, dec: float, tolerance_deg=0.01) -> StarDescription:
        logging.debug(f"get_ucac4_sd with ra:{ra}, dec:{dec}, tolerance:{tolerance_deg}")
        # don't want to bother with more than 2 zone overlappings
        target_np = np.array((ra, dec))
        assert tolerance_deg < 0.2

        zone1 = self.get_zone_for_dec(dec - tolerance_deg / 2)
        zone2 = self.get_zone_for_dec(dec + tolerance_deg / 2)
        # check for zone overlap
        zones = [zone1]
        if zone1 != zone2:
            zones.append(zone2)
        logging.debug(f"Zone range is {zones}")
        bucket1 = self.get_ra_bucket(ra - tolerance_deg / 2)
        bucket2 = self.get_ra_bucket(ra + tolerance_deg / 2)
        buckets = range(bucket1, bucket2 + 1)
        logging.debug(f"Bucket range is from {bucket1} to {bucket2}")
        smallest_dist = 1000
        best = None
        for zone in zones:
            for bucket in buckets:
                cache_key = (zone, bucket)
                cache = self.zone_bucket_cache.get(cache_key)
                if cache != -1:
                    the_stars = cache
                else:
                    logging.debug(f"Processing {zone}:{bucket}")
                    index = (bucket - 1) * 900 + zone - 1
                    star_run_nr = self.result_n0[index]
                    star_count_in_bucket = self.result_nn[index]
                    the_range = list(range(star_run_nr, star_run_nr + star_count_in_bucket))
                    the_stars = self.get_ucac4_details_raw(zone, the_range)
                    self.zone_bucket_cache.set(cache_key, the_stars)
                for sd in the_stars:
                    dist = np.linalg.norm(np.array((UCAC4.get_real_ra_dec(sd[0].ra, sd[0].spd))-target_np))
                    if dist < smallest_dist:
                        smallest_dist = dist
                        best = sd
        return self.get_ucac4_star_description_fromtuple(*best)


    @staticmethod
    def get_zone_for_dec(dec: float):
        dec_0 = dec + 90
        index = round(dec_0 / 0.2, 1)
        return math.ceil(index) if index != 0 else 1


    @staticmethod
    def get_ra_bucket(ra: float):
        index = math.ceil(ra * 4)  # 1440/360
        return index if index != 0 else 1


    @staticmethod
    def name_to_zone_and_run_nr(ucac_id: str):
        full_id = ucac_id[6:]
        zone = full_id[:3]
        run_nr = int(full_id[4:].lstrip("0"))
        return zone, run_nr


    @staticmethod
    def zone_and_run_nr_to_name(zone: str, run_nr: int):
        return f"UCAC4 {zone}-{run_nr:06}"


    def get_ucac4_star_description_fromid(self, ucac4_id) -> StarDescription:
        return UCAC4.get_ucac4_star_descriptions(self.get_ucac4_details(ucac4_id))[0]


    @staticmethod
    def get_ucac4_star_descriptions(startuples: List[Tuple[StarTuple, str, int]]) -> List[StarDescription]:
        result = []
        for startuple in startuples:
            star, zone, run_nr = startuple
            result.append(UCAC4.get_ucac4_star_description_fromtuple(star, zone, run_nr))
        return result


    @staticmethod
    def get_ucac4_star_description_fromtuple(star: StarTuple, zone: str, run_nr: int):
        ra, dec = UCAC4.get_real_ra_dec(star.ra, star.spd)
        sd = StarDescription(coords=SkyCoord(ra, dec, unit='deg'),
                             vmag=star.apass_mag_V / 1000, e_vmag=abs(star.apass_mag_sigma_V / 100),
                             aavso_id=UCAC4.zone_and_run_nr_to_name(zone, run_nr))
        return sd

    @staticmethod
    def get_real_ra_dec(ra, spd) -> Tuple[float, float]:
        return ra / 1000 / 3600, (spd - 324000000) / 1000 / 3600

    @staticmethod
    def make_star(result: List):
        star = StarTuple._make(result)
        star = star._replace(ra_sigma=star.ra_sigma + 128)
        star = star._replace(dec_sigma=star.dec_sigma + 128)
        star = star._replace(pm_ra_sigma=star.pm_ra_sigma + 128)
        star = star._replace(pm_dec_sigma=star.pm_dec_sigma + 128)
        return star

    def add_ucac4_to_sd(self, stars: List[StarDescription]):
        with tqdm.tqdm(total=len(stars), desc='Adding UCAC4', unit='stars') as pbar:
            for star in stars:
                if not star.has_metadata("UCAC4"):
                    sd = self.get_ucac4_sd(star.coords.ra.deg, star.coords.dec.deg)
                    do_calibration.add_info_to_star_description(star, sd.vmag, sd.e_vmag, sd.aavso_id, "UCAC4",
                                                                sd.coords)
                pbar.update(1)

    # >>> ra=140361
    # >>> ra/1000/60/60
    # 0.038989166666666665
    # >>> dec=324005500
    # >>> (dec-324000000)/1000/60/60
    # 0.0015277777777777776
    # >>>




if __name__ == '__main__':
    # run('UCAC4 003-000419')
    ucac4 = UCAC4()
    print(ucac4.get_ucac4_star_description('UCAC4 001-000003'))

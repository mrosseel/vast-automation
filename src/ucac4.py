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
import os
import logging
import math
import struct
import numpy as np
from collections import namedtuple
from typing import List, Tuple, Dict
from star_description import StarDescription
from astropy.coordinates import SkyCoord
from pathlib import Path
from LRUCache import LRUCache

# remove this, would be cleaner in ucac4_utils or something
import do_calibration
import tqdm
import decimal

StarTuple = namedtuple(
    "Star",
    "ra spd mag1 mag2 mag_sigma obj_type double_star_flag ra_sigma"
    " dec_sigma n_ucac_total n_ucac_used n_cats_used epoch_ra epoch_dec pm_ra pm_dec"
    " pm_ra_sigma pm_dec_sigma twomass_id mag_j mag_h mag_k icq_flag1 icq_flag2 icq_flag3"
    " e2mpho1 e2mpho2 e2mpho3 apass_mag_B apass_mag_V apass_mag_g apass_mag_r apass_mag_i"
    " apass_mag_sigma_B apass_mag_sigma_V apass_mag_sigma_g apass_mag_sigma_r"
    " apass_mag_sigma_i yale_gc_flags catalog_flags leda_flag twomass_ext_flag"
    " id_number ucac2_zone ucac2_number",
)
MinimalStarTuple = namedtuple("MinStar", "id, ra, dec, mag")
UcacTuple = Tuple[StarTuple, str, int]
UcacTupleList = List[UcacTuple]


def get_line_nr(n0, nn, line):
    return f"{n0[line * 900]}\t{nn[line * 900]}"


# 1  n0   = running star number (index along the main data file) of the star before the first one in this bin,
#           the sequence starts out with 0 at the beginning of each new declination zone
# 2  nn   = number of stars in this bin (which can be zero)
# 3  zn   = zone number (1 to 900)
# 4  j    = index for bins along RA (1 to 1440)
# 5  dec  = upper declination of corresponding zone, printed out only at the beginning of a new zone


class UCAC4:
    """ Uses UCAC4 catalog to search ucac4 numbers from coords and vice versa. Ucac4 path can be passed or read from $ucac4_path """

    def __init__(self, ucac_path=None):
        if ucac_path is None:
            try:
                ucac_path = Path(os.environ["ucac4_path"])
                print("found environ with ucac45 path", ucac_path)
            except KeyError:
                logging.debug(
                    "No ucac path passed, and no environment var, using default"
                )
                ucac_path = Path("./support/ucac4/UCAC4/")
        logging.info("Ucac path is", ucac_path)
        # star id given by munipack
        self.ucac_path = ucac_path
        self.zones_cache = LRUCache(capacity=10)
        self.bucket_cache = LRUCache(capacity=100000)
        self.zone_bucket_cache = LRUCache(capacity=10000)
        self.index_cache = None
        self.ra_range = 360 * 3600 * 100
        # read index file
        with open(
            str(Path(ucac_path, "u4i", "u4index.unf")), mode="rb"
        ) as file:  # b is important -> binary
            self.index_cache = file.read()
        self.zone_starformat = "=iiHHBBBbbBBBHHhhbbIHHHBBBBBBHHHHHbbbbbBIBBIHI"
        self.zone_star_length = struct.calcsize(self.zone_starformat)
        self.unpack_zone_fileformat = struct.Struct(self.zone_starformat).unpack
        (
            self.result_n0_running_star_number,
            self.result_nn_stars_in_bin,
        ) = UCAC4._get_n0_and_nn(self.index_cache)

    def get_zone_filecontent(self, zone: int):
        """ gets the content of a zone file, either from disk or from cache"""
        result = self.zones_cache.get(zone)
        if result != -1:
            return result
        with open(
            str(Path(self.ucac_path, f"u4b/z{zone:03}")), mode="rb"
        ) as file:  # b is important -> binary
            result = file.read()
            self.zones_cache.set(zone, result)
            return result

    def get_ucactuple_from_id(self, ucac_id) -> UcacTuple:
        """ Given a UCAC ID, return a tuple of (StarTuple, zone, run_nr) """
        zone, run_nr = UCAC4.ucac_id_to_zone_and_run_nr(ucac_id)
        logging.debug(f"UCAC4 id {zone}, {run_nr}")
        return self.get_ucactuples_for_zone_and_runnrs(zone, [run_nr])[0]

    def index_bin_to_run_nrs(self, zone: int, index_bin: int):
        # index = (zone - 1) * 1440 + index_bin
        index = (index_bin - 1) * 900 + zone - 1
        star_run_nr = self.result_n0_running_star_number[index]
        star_count_in_bucket = self.result_nn_stars_in_bin[index]
        run_nrs = list(range(star_run_nr, star_run_nr + star_count_in_bucket))
        logging.debug(
            f"index_bin_to_run_nrs: zone is {zone}, index_bin = {index_bin}, index is {index}. "
            f"star_run_nr is {star_run_nr}, count =  {star_count_in_bucket}, run nrs: {run_nrs}"
        )
        return run_nrs

    def get_zones_and_index_bins(self, ra, dec, tolerance_deg) -> Dict[int, List[int]]:
        logging.debug(f"ra: {ra}, dec: {dec}, tolerance: {tolerance_deg}")
        zone = self.get_zone_for_dec(dec - tolerance_deg / 2)
        end_zone = self.get_zone_for_dec(dec + tolerance_deg / 2)
        # check for zone overlap
        zones = list(range(zone, end_zone + 1))
        logging.debug(f"get_zones_and_index_bins: Zone range is {zones}")
        ra_low = self.ra_bin_index(ra - tolerance_deg / 2)  # ra_start in C code
        ra_high = self.ra_bin_index(ra + tolerance_deg / 2)
        index_bins = list(range(ra_low, ra_high + 1))
        logging.debug(f"get_zones_and_index_bins: {ra_low}, {ra_high}, {index_bins}")
        return zones, index_bins

    def get_ucactuples_for_zone_and_runnrs(
        self, zone: int, run_nr_list: List[int]
    ) -> UcacTupleList:
        """ Given a zone and one or more run_nr's, return List of (StarTuple, zone, run_nr) """
        stars = []
        filecontent = self.get_zone_filecontent(zone)
        for run_nr in run_nr_list:
            key = (zone, run_nr)
            # logging.debug(f"Getting star from zone {zone}, run_nr {run_nr}")
            star = self.bucket_cache.get(key)
            if star == -1:
                result = self.unpack_zone_fileformat(
                    filecontent[
                        self.zone_star_length
                        * (run_nr - 1) : self.zone_star_length
                        * run_nr
                    ]
                )
                star = UCAC4.make_startuple(result)
                self.bucket_cache.set(key, star)
            stars.append((star, zone, run_nr))
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

    def get_region_minimal_star_tuples(
        self, ra: float, dec: float, radius=0.5
    ) -> List[MinimalStarTuple]:
        """ For a given ra/dec and radius, return all ucac4 stars as ('id, ra, dec, mag') """
        zones, buckets = self.get_zones_and_index_bins(ra, dec, radius)
        result = []
        for zone in zones:
            for bucket in buckets:
                cache_key = (zone, bucket)
                ucactuples: UcacTupleList = self.zone_bucket_cache.get(cache_key)
                # cache miss
                if ucactuples == -1:
                    ucactuples: UcacTupleList = self.get_ucactuples_for_zone_and_runnrs(
                        zone, self.index_bin_to_run_nrs(zone, bucket)
                    )
                    self.zone_bucket_cache.set(cache_key, ucactuples)
                if len(ucactuples) == 0:
                    logging.debug(f"zone/bucket: {zone}/{bucket}, no stars")
                    if bucket + 1 not in buckets:
                        buckets.append(bucket + 1)
                        logging.debug(f"Appending bucket {bucket+1}")
                for ucactuple in ucactuples:
                    result.append(
                        MinimalStarTuple(
                            UCAC4.zone_and_run_nr_to_name(ucactuple[1], ucactuple[2]),
                            *UCAC4.get_real_ra_dec(ucactuple[0].ra, ucactuple[0].spd),
                            ucactuple[0].apass_mag_V / 1000,
                        )
                    )
        return result

    def get_sd_from_ra_dec(
        self, ra: float, dec: float, tolerance_deg=0.02
    ) -> StarDescription:
        return self.get_star_description_from_tuple(
            self.get_ucactuple_from_ra_dec(ra, dec, tolerance_deg)
        )

    def get_ucactuple_from_ra_dec(
        self, ra: float, dec: float, tolerance_deg=0.02
    ) -> UcacTuple:
        logging.debug(
            f"get_ucac4_ucactuple_from_ra_dec with ra:{ra}, dec:{dec}, tolerance:{tolerance_deg}"
        )
        target_np = np.array((ra, dec))
        zones, buckets = self.get_zones_and_index_bins(ra, dec, tolerance_deg)
        smallest_dist = 1000
        best = None
        for zone in zones:
            for bucket in buckets:
                cache_key = (zone, bucket)
                ucactuples: UcacTupleList = self.zone_bucket_cache.get(cache_key)
                # cache miss
                if ucactuples == -1:
                    ucactuples: UcacTupleList = self.get_ucactuples_for_zone_and_runnrs(
                        zone, self.index_bin_to_run_nrs(zone, bucket)
                    )
                    self.zone_bucket_cache.set(cache_key, ucactuples)
                ################ DEBUG
                if len(ucactuples) > 0:
                    radecs = [
                        self.get_real_ra_dec(x[0].ra, x[0].spd) for x in ucactuples
                    ]
                    ras = [x[0] for x in radecs]
                    decs = [x[1] for x in radecs]
                    logging.debug(
                        f"zone/bucket: {zone}/{bucket}, Searching between {min(ras)}, {max(ras)}, {min(decs)}, {max(decs)}"
                    )
                else:
                    logging.debug(f"zone/bucket: {zone}/{bucket}, no stars")
                    if bucket + 1 not in buckets:
                        buckets.append(bucket + 1)
                        logging.debug(f"Appending bucket {bucket+1}")
                ################ DEBUG
                for ucactuple in ucactuples:
                    dist = np.linalg.norm(
                        np.array(
                            (UCAC4.get_real_ra_dec(ucactuple[0].ra, ucactuple[0].spd))
                            - target_np
                        )
                    )
                    # logging.debug(f"magj: {sd[0].mag_j}")
                    if dist < smallest_dist:
                        smallest_dist = dist
                        best = ucactuple
        if best is None:
            logging.warning(
                f"Did not find a UCAC4 match for {ra}, {dec}, {tolerance_deg}. Buckets: {buckets}, "
                f"zones: {zones},smallest dist: {smallest_dist}"
            )
            return best
        logging.debug(f"Best distance is: {smallest_dist}, {best}")
        return best

    @staticmethod
    def get_zone_for_dec(dec: float, zone_height: float = 0.2) -> int:
        dec_0 = decimal.Decimal(dec) + decimal.Decimal(90.0)
        result = min(900, max(1, int(dec_0 / decimal.Decimal(zone_height)) + 1))
        logging.debug(
            f"get_zone_for_dec: dec {dec}, height:{zone_height}, dec0 {dec_0}, result {result}"
        )
        return result

    @staticmethod
    def ra_bin_index(ra: float):
        """index for bins along RA (1 to 1440)"""
        index = math.ceil(decimal.Decimal(ra) * decimal.Decimal(4))  # 1440/360
        logging.debug(
            f"ra_bin_index: ra {ra}, index {index}, rawindex: {decimal.Decimal(ra) * decimal.Decimal(4)}"
        )
        return max(1, index)

    @staticmethod
    def ucac_id_to_zone_and_run_nr(ucac_id: str):
        full_id = ucac_id[6:]
        zone = int(full_id[:3])
        run_nr = int(full_id[4:].lstrip("0"))
        return zone, run_nr

    @staticmethod
    def zone_and_run_nr_to_name(zone: int, run_nr: int):
        return f"UCAC4 {zone:03}-{run_nr:06}"

    def get_star_description_from_id(self, ucac4_id) -> StarDescription:
        return UCAC4.get_star_description_from_tuple(
            self.get_ucactuple_from_id(ucac4_id)
        )

    @staticmethod
    def get_star_descriptions_from_tuples(
        ucactuples: UcacTupleList,
    ) -> List[StarDescription]:
        return [
            UCAC4.get_star_description_from_tuple(ucactuple) for ucactuple in ucactuples
        ]

    @staticmethod
    def get_star_description_from_tuple(ucactuple: UcacTuple) -> StarDescription:
        startuple, zone, run_nr = ucactuple
        ra, dec = UCAC4.get_real_ra_dec(startuple.ra, startuple.spd)
        sd = StarDescription(
            coords=SkyCoord(ra, dec, unit="deg"),
            vmag=startuple.apass_mag_V / 1000,
            vmag_err=abs(startuple.apass_mag_sigma_V / 100),
            aavso_id=UCAC4.zone_and_run_nr_to_name(zone, run_nr),
        )
        return sd

    @staticmethod
    def get_real_ra_dec(ra, spd) -> Tuple[float, float]:
        # return ra / 1000 / 3600, (spd - 324000000) / 1000 / 3600
        divisor = decimal.Decimal(3600000)
        return (
            float(decimal.Decimal(ra) / divisor),
            float(decimal.Decimal(spd - 324000000) / divisor),
        )

    @staticmethod
    def make_startuple(result: List) -> StarTuple:
        star = StarTuple._make(result)
        star = star._replace(ra_sigma=star.ra_sigma + 128)
        star = star._replace(dec_sigma=star.dec_sigma + 128)
        star = star._replace(pm_ra_sigma=star.pm_ra_sigma + 128)
        star = star._replace(pm_dec_sigma=star.pm_dec_sigma + 128)
        return star

    def add_sd_metadatas(self, stars: List[StarDescription], overwrite=False):
        with tqdm.tqdm(total=len(stars), desc="Adding UCAC4", unit="stars") as pbar:
            for star in stars:
                if not star.has_metadata("UCAC4") or overwrite:
                    sd = self.get_sd_from_ra_dec(
                        star.coords.ra.deg, star.coords.dec.deg
                    )
                    self._add_catalog_data_to_sd(star, sd, overwrite)
                pbar.update(1)

    def add_sd_metadata_from_id(
        self, star: StarDescription, ucac4_id: str, overwrite=False
    ):
        sd = self.get_star_description_from_id(ucac4_id)
        self._add_catalog_data_to_sd(star, sd, overwrite)

    @staticmethod
    def _add_catalog_data_to_sd(
        sd: StarDescription, ucac4_sd: StarDescription, overwrite
    ):
        """ Add UCAC4 catalog data to a stardescription if there is none yet, or if overwrite is True """
        if ucac4_sd is not None and not sd.has_metadata("UCAC4") or overwrite:
            do_calibration.add_catalog_data_to_sd(
                sd,
                ucac4_sd.vmag,
                ucac4_sd.vmag_err,
                ucac4_sd.aavso_id,
                "UCAC4",
                ucac4_sd.coords,
            )

    # >>> ra=140361
    # >>> ra/1000/60/60
    # 0.038989166666666665
    # >>> dec=324005500
    # >>> (dec-324000000)/1000/60/60
    # 0.0015277777777777776
    # >>>


if __name__ == "__main__":
    # run('UCAC4 003-000419')
    ucac4 = UCAC4()
    print(ucac4.get_ucac4_star_description("UCAC4 001-000003"))

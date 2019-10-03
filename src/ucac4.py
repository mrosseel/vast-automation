
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
import struct
from collections import namedtuple
from typing import List
from star_description import StarDescription
from astropy.coordinates import SkyCoord


Star = namedtuple('Star', 'ra spd mag1 mag2 mag_sigma obj_type double_star_flag ra_sigma' 
                          ' dec_sigma n_ucac_total n_ucac_used n_cats_used epoch_ra epoch_dec pm_ra pm_dec' 
                          ' pm_ra_sigma pm_dec_sigma twomass_id mag_j mag_h mag_k icq_flag1 icq_flag2 icq_flag3'
                          ' e2mpho1 e2mpho2 e2mpho3 apass_mag_B apass_mag_V apass_mag_g apass_mag_r apass_mag_i'
                          ' apass_mag_sigma_B apass_mag_sigma_V apass_mag_sigma_g apass_mag_sigma_r apass_mag_sigma_i'
                          ' yale_gc_flags catalog_flags leda_flag twomass_ext_flag id_number ucac2_zone ucac2_number')

def get_ucac4_details(ucac_id):
    id = ucac_id[6:]
    zone = id[:3]
    run_nr = int(id[4:].lstrip("0"))

    print("id", id, zone, run_nr)
    with open(f"./support/ucac4/UCAC4/u4b/z{zone}" , mode='rb') as file:  # b is important -> binary
        fileContent = file.read()
        print("filecontent mod 78 is ", len(fileContent) % 78)
        # print("filecontent", fileContent[:78])
        # starformat = "=iihhBBBBBBBBhhhhBBihhhBBBBBBhhhhhBBBBBBiBBihi"
        starformat = "=iiHHBBBbbBBBHHhhbbIHHHBBBBBBHHHHHbbbbbBIBBIHI"
        starlength = struct.calcsize(starformat)
        print("starlength is", starlength)
        unpack_starformat = struct.Struct(starformat).unpack
        result = unpack_starformat(fileContent[starlength*(run_nr-1):starlength*run_nr])
        star = make_star(result)
        print("read ucac4 star:", star)
        return star

def get_ucac4_star_description(ucac4_id):
    star = get_ucac4_details(ucac4_id)
    sd = StarDescription(coords=SkyCoord(star.ra/1000/3600, (star.spd-324000000)/1000/3600, unit='deg'),
                           vmag=star.apass_mag_V/1000, e_vmag=star.apass_mag_sigma_V/100, aavso_id=ucac4_id)
    return sd


def make_star(result: List):
    star = Star._make(result)
    star = star._replace(ra_sigma=star.ra_sigma+128)
    star = star._replace(dec_sigma=star.dec_sigma+128)
    star = star._replace(pm_ra_sigma=star.pm_ra_sigma+128)
    star = star._replace(pm_dec_sigma=star.pm_dec_sigma+128)
    return star

# >>> ra=140361
# >>> ra/1000/60/60
# 0.038989166666666665
# >>> dec=324005500
# >>> (dec-324000000)/1000/60/60
# 0.0015277777777777776
# >>>

if __name__ == '__main__':
    # run('UCAC4 003-000419')
    print(get_ucac4_star_description('UCAC4 001-000003'))

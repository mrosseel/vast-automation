from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pandas as pd
import os
import init
import reading

def calculate_wcs_manual(reference_frame, xpos, ypos, arcsec_width, arcsec_height):
    object_ra, object_dec, naxis1, naxis2, jd= getDataFromFitsHeader(reference_frame)
    c = SkyCoord(object_ra, object_dec, unit="deg")
    return [setup_wcs(c, naxis1, naxis2, xpos, ypos, arcsec_width, arcsec_height), jd]

def calculate_wcs_from_file(header, frame, xpos, ypos):
    object_ra, object_dec, naxis1, naxis2, jd= getDataFromFitsHeader(frame)
    c = SkyCoord(object_ra, object_dec, unit="deg")
    return [setup_wcs_from_file(header, c, xpos, ypos), jd]


def setup_wcs(coord, naxis1, naxis2, xpos, ypos, arcsec_width, arcsec_height):
    w = wcs.WCS(naxis=2)
    # Set up an "Airy's zenithal" projection
    # Vector properties may be set with Python lists, or Numpy arrays
    w.wcs.crpix = [xpos, ypos] # the so-called center, this is the position of our main star, given by the user
    w.wcs.cdelt = np.array([-0.000572222222222, 0.000572222222222])
    w.wcs.crval = [coord.ra.degree, coord.dec.degree]
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    return w

def setup_wcs_from_file(header, coord, xpos, ypos):
    header_fits = get_fits_header(header)
    w = wcs.WCS(header_fits)
    w.wcs.crpix = [xpos, ypos] # the so-called center, this is the position of our main star, given by the user
    w.wcs.crval = [coord.ra.degree, coord.dec.degree]
    return w

### HELPER functions

def getDataFromFitsHeader(reference_frame):
    fits_header = get_fits_header(reference_frame)
    object_ra = fits_header['OBJCTRA']
    object_dec = fits_header['OBJCTDEC']
    naxis1 = fits_header['NAXIS1']
    naxis2 = fits_header['NAXIS2']
    jd = fits_header['JD']
    return [object_ra, object_dec, naxis1, naxis2, jd]

def pixel_to_radec(wcs_config, xpix, ypix):
    pixcrd = np.array([[xpix, ypix]], np.float_)
    result = wcs_config.wcs_pix2world(pixcrd, 1)
    return SkyCoord(result[0][0], result[0][1], unit='deg')

def star_to_radec(star, w, jd):
    jd, x, y, mag = reading.read_pos(star, jd)
    print("star to radec head", jd, x, y, mag)
    radec = pixel_to_radec(w, x, y)
    return [radec, mag]

def get_fits_header(reference_file):
    hdulist = fits.open(reference_file)
    return hdulist[0].header

### TEST code

def testing(w):
    # Some pixel coordinates of interest.
    pixcrd = np.array([[0, 0], [1365/2.0, 1365/2.0], [0, 1365]], np.float_)
    print(pixcrd)

    # Convert pixel coordinates to world coordinates
    world = w.wcs_pix2world(pixcrd, 1)
    print(world)

    # Convert the same coordinates back to pixel coordinates.
    pixcrd2 = w.wcs_world2pix(world, 1)
    print(pixcrd2)

    # These should be the same as the original pixel coordinates, modulo
    # some floating-point error.
    assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6
    print("assert ok")


def test_code(reference_frame, xpos, ypos, arcsecond_width, arssecond_heigth):
    reference_frame_full = init.basedir +reference_frame
    fits_header = get_fits_header(reference_frame_full)
    object_ra = fits_header['OBJCTRA']
    object_dec = fits_header['OBJCTDEC']
    naxis1 = fits_header['NAXIS1']
    naxis2 = fits_header['NAXIS2']
    JD = fits_header['JD']

    print("object_ra:", object_ra, "object_dec:", object_dec, "naxis1:", naxis1, "naxis2:", naxis2)

    # coords of center of frame of first image
    c = SkyCoord(object_ra, object_dec, unit="deg")
    #wcs_config = setup_wcs(c, naxis1, naxis2)
    #print("file:",reference_frame, os.getcwd())
    #wcs_config = setup_wcs_from_file(init.basedir +reference_frame)
    #result = pixel_to_radec(wcs_config, 1365, 1365)
    #print(result)

    # conesearch is deprecated and moved to astroquery
    #conesearch.list_catalogs()
    #my_catalog = 'Guide Star Catalog v2 1'
    #c = SkyCoord.from_name('GSC 7911-3668')
    #result = conesearch.conesearch(c, 0.01 * u.degree, catalog_db=my_catalog)
    #print(result)


    # Set the WCS information manually by setting properties of the WCS
    # object.

    # Create a new WCS object.  The number of axes must be set
    # from the start
    w = wcs.WCS(naxis=2)
    print(c )

    # Set up an "Airy's zenithal" projection
    # Vector properties may be set with Python lists, or Numpy arrays
    w.wcs.crpix = [xpos, ypos] # the so-called center, this is the position of our main star, given by the user
    w.wcs.cdelt = np.array([-0.000572222222222, 0.000572222222222])
    #w.wcs.crval = [spaces_position_to_comma(object_ra), spaces_position_to_comma(object_dec)]
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]

    # Some pixel coordinates of interest.
    pixcrd = np.array([[0, 0], [1365/2.0, 1365/2.0], [0, 1365]], np.float_)
    print(pixcrd)

    # Convert pixel coordinates to world coordinates
    world = w.wcs_pix2world(pixcrd, 1)
    print(world)

    # Convert the same coordinates back to pixel coordinates.
    pixcrd2 = w.wcs_world2pix(world, 1)
    print(pixcrd2)

    # These should be the same as the original pixel coordinates, modulo
    # some floating-point error.
    assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6
    print("assert ok")

    # Now, write out the WCS object as a FITS header
    header = w.to_header()

    # header is an astropy.io.fits.Header object.  We can use it to create a new
    # PrimaryHDU and write it to a file.
    hdu = fits.PrimaryHDU(header=header)
    # Save to FITS file
    # hdu.writeto('test.fits')

#w, jd = calculate_wcs('WWCrA#30V_000185980_FLAT.fit', 695, 671, 47*60, 47*60)
w, jd = calculate_wcs_from_file(init.reference_header, init.reference_frame, init.xpos, init.ypos)
print(w, jd)
testing(w)
pd.set_option('precision', 10)
print(star_to_radec(1, w, jd))

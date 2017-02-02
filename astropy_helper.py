import init # should go
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.vo.client import conesearch
from astropy import units as u

#cell 5
def setup_wcs(coord, naxis1, naxis2):
    w = wcs.WCS(naxis=2)

    # Set up an "Airy's zenithal" projection
    # Vector properties may be set with Python lists, or Numpy arrays
    w.wcs.crpix = [naxis1/2.0, naxis2/2.0]
    w.wcs.cdelt = np.array([-0.000572222222222, 0.000572222222222])
    print("ra/dec", coord.ra.degree, coord.dec.degree)
    w.wcs.crval = [coord.ra.degree, coord.dec.degree]
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    return w

def setup_wcs_from_file(filename):
    hdulist = fits.open(filename)
    return wcs.WCS(hdulist[0].header)

def pixel_to_radec(wcs_config, xpix, ypix):
    pixcrd = np.array([[xpix, ypix]], np.float_)
    result = wcs_config.wcs_pix2world(pixcrd, 1)
    return SkyCoord(result[0][0], result[0][1], unit='deg')

def star_to_radec():
    print("todo")


def get_fits_header(reference_file):
    hdulist = fits.open(reference_file)
    return hdulist[0].header


def dit_was_in_de_root():
    # fits_header = get_fits_header(reference_frame)
    # object_ra = fits_header['OBJCTRA']
    # object_dec = fits_header['OBJCTDEC']
    # naxis1 = fits_header['NAXIS1']
    # naxis2 = fits_header['NAXIS2']


    # coords of center of frame of first image
    # c = SkyCoord(object_ra, object_dec, unit="deg")
    #wcs_config = setup_wcs(c, naxis1, naxis2)
    print("file:",reference_frame, os.getcwd())
    wcs_config = setup_wcs_from_file('./'+reference_frame)
    result = pixel_to_radec(wcs_config, 1365, 1365)
    print(result)

    #cell 6
    #conesearch.list_catalogs()
    my_catalog = 'Guide Star Catalog v2 1'
    #c = SkyCoord.from_name('GSC 7911-3668')
    #result = conesearch.conesearch(c, 0.01 * u.degree, catalog_db=my_catalog)
    print(result)


    #cell 7
    # Set the WCS information manually by setting properties of the WCS
    # object.

    # Create a new WCS object.  The number of axes must be set
    # from the start
    w = wcs.WCS(naxis=2)

    # Set up an "Airy's zenithal" projection
    # Vector properties may be set with Python lists, or Numpy arrays
    w.wcs.crpix = [1365/2.0, 1365/2.0]
    w.wcs.cdelt = np.array([-0.000572222222222, 0.000572222222222])
    w.wcs.crval = [18.0936111, -43.8325000]
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]

    # Some pixel coordinates of interest.
    pixcrd = np.array([[0, 0], [1365/2.0, 1365/2.0], [45, 98]], np.float_)

    # Convert pixel coordinates to world coordinates
    world = w.wcs_pix2world(pixcrd, 1)
    print(world)

    # Convert the same coordinates back to pixel coordinates.
    pixcrd2 = w.wcs_world2pix(world, 1)
    print(pixcrd2)

    # These should be the same as the original pixel coordinates, modulo
    # some floating-point error.
    assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6

    # Now, write out the WCS object as a FITS header
    header = w.to_header()

    # header is an astropy.io.fits.Header object.  We can use it to create a new
    # PrimaryHDU and write it to a file.
    hdu = fits.PrimaryHDU(header=header)
    # Save to FITS file
    # hdu.writeto('test.fits')

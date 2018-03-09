import init
from astropy.wcs import WCS

def calibrate():
    w = WCS(init.basedir + 'new-image.fits')
    #xpos, ypos = w.all_world2pix(init.ra_deg, init.dec_deg, 0, ra_dec_order=True)
    #print(xpos, ypos)
    #print(w)
    return w

calibrate()
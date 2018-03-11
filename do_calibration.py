import os
import init
from astropy.wcs import WCS

def calibrate():
    w = WCS(init.reference_header)
    #xpos, ypos = w.all_world2pix(init.ra_deg, init.dec_deg, 0, ra_dec_order=True)
    #print(xpos, ypos)
    #print(w)
    return w

def find_reference_in_files(the_path):
    the_dir = os.listdir(the_path)
    the_dir.sort()
    reference_frame_index = the_dir.index(init.reference_frame)
    print(reference_frame_index)

find_reference_in_files(init.fitsdir)
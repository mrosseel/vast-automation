# VAST doesn't like strange chars like '#'. This file will delete them.
import os
import sys


def doit(thepath):
    path = thepath 
    print(f"Renaming all files in {path}")
    for filename in os.listdir(path):
        newfilename = filename.replace('#', '_')
        print(f"Renaming {filename} to {newfilename}")
        os.rename(os.path.join(path,filename), os.path.join(path,newfilename))

if __name__ == '__main__':
    doit(sys.argv[1])

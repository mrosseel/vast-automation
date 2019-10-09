# munipack-automation

Docker is used to construct a light-weight virtual machine linux containing all necessary dependencies.
Once you're in this VM, all python commands can be run.

## Preparation

### Importing VSX star catalog

* AAVSO VSX catalog can be downloaded here: http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=B%2Fvsx%2Fversions%2F2018-02-26&target=brief&
* run 'python vsx_pickle.py vsx.dat' where vsx.dat is the unzipped versin of the downloaded vsx catalog
* check that 'vsx_catalog.bin' has been written successfully

### Starting Docker

* have a working docker installation: https://www.docker.com/community-edition
* cd docker
* docker build . -t mrosseel/munipack-automation
* cd ..
* command: `./startJupyter.sh`
* you are automatically logged into a root shell of the docker container

## Run VAST on the fits files

`./vast -u -x 3 ../location/of/fits/*.fit`

This uses UTC time, and ignores the 'blended' flag for stars which are close to each other.
This will generate many vast files in the vast directory

## Process VAST results

### run the command

`./vast_process.sh --vsx --candidates -d support/vast-1.0rc84`

### plate solve the reference frame (first run only)

The software will stop and ask you to do this.

* take the reference frame and calculate a fits header using http://Astrometry.net
* store it in the vast directory as *new-image.fits*

### look at the results

This command line above will generate vsx information and create phase diagrams for 
all vast autocandidates and vsx stars.
Also a few extra files are generated:

* vsx_stars.txt
* vast_list_of_all_stars_pos.txt
* vast_autocandidates_pos.txt

## Other

### Jupyter lab usage

The docker image also exposes a Jupyter lab instance on port 8888.
_Password is 'muni'_

### Importing UCAC4 star catalog

Getting the 900 files (9Gb):
- `wget ftp://cdsarc.u-strasbg.fr/0/more/UCAC4/u4b/*`

Checking that all 900 files were downloaded correctly:
- `md5sum -c md5sum.txt`

## TODO

- calculate correct realMag, realErr
- for candidates, closest known vsx
- for candidates_pos, no h/m/s but spaces (see get_lesve)
+ check that if no comparison stars are found 0 is used
- make html page with selected candidates
/ make lightcurves based on JD
+ hire kimsufi server for processing ks-11
- file to give comparison stars -f bla.txt
- check out https://public.lanl.gov/palmer/fastchi.html for period determination

- AAVSO report should use instrumental magnitudes for comparison stars
- other stars in command line
- check out https://github.com/toros-astro/astroalign for aligning
- new ensemble comparison star calculation?
- stacking images to detect fainter stars+have better signal/noise ratio: https://github.com/fedhere/coaddfitim

## References

* https://books.google.ae/books?id=g_K3-bQ8lTUC&pg=PA236&lpg=PA236&dq=munipack+aperture&source=bl&ots=P4BKKI25HG&sig=Lj9Kg6EZi2pwKZXK5Hk5_B4qcVg&hl=en&sa=X&ved=0ahUKEwj4qr6o4evZAhWILcAKHdtnAqYQ6AEIRDAE#v=onepage&q&f=false
* https://www.aavso.org/how-report-new-variable-star-discoveries

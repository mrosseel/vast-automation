# munipack-automation

Docker is used to construct a light-weight virtual machine linux containing all necessary dependencies.
Once you're in this VM, all python commands can be run.

## Starting Docker

* have a working docker installation: https://www.docker.com/community-edition
* cd docker
* docker build . -t mrosseel/munipack-automation
* cd ..
* command: `./startJupyter.sh`
* you are automatically logged into a root shell of the docker container

## Command line usage

### Init settings

* take one reference frame and calculate a fits header using http://Astrometry.net
* edit init.py to set all correct directories and values

### Run

* command: `python do_muniwin.py`

### File overview

* init.py : directory settings, processing settings
* do_muniwin.py : start all
* do_charts.py : plots of lightcurve and phase diagrams
* do_field_charts.py : plot the reference frame + circles around stars of interest
* do_upsilon.py : only do machine learning detection
* do_aavso_report: write files in the AAVSO Extended Format
* do_profile.py : do performance profiling on the app (not sure if working)

## Jupyter lab usage

The docker image also exposes a Jupyter lab instance on port 8888.
_Password is 'muni'_

## Importing VSX star catalog

* AAVSO VSX catalog can be downloaded here: http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=B%2Fvsx%2Fversions%2F2018-02-26&target=brief&
* run 'python vsx_pickle.py vsx.dat' where vsx.dat is the unzipped versin of the downloaded vsx catalog
* check that 'vsx_catalog.bin' has been written succesfully

## TODO

- check out https://public.lanl.gov/palmer/fastchi.html for period determination
- new ensemble comparison star calculation?
- check error column for any stars having error bars > 1%
- stacking images to detect fainter stars+have better signal/noise ratio: https://github.com/fedhere/coaddfitim
- write munifind_ intermediate files in a folder
- detect less stars to do the aperture calculations (generate them in a seperate folder, subset of images)

## References

* https://books.google.ae/books?id=g_K3-bQ8lTUC&pg=PA236&lpg=PA236&dq=munipack+aperture&source=bl&ots=P4BKKI25HG&sig=Lj9Kg6EZi2pwKZXK5Hk5_B4qcVg&hl=en&sa=X&ved=0ahUKEwj4qr6o4evZAhWILcAKHdtnAqYQ6AEIRDAE#v=onepage&q&f=false
* https://www.aavso.org/how-report-new-variable-star-discoveries

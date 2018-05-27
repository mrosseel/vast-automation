# munipack-automation

Docker is used to construct a light-weight virtual machine linux containing all necessary dependencies.
Once you're in this VM, all python commands can be run.

## Starting Docker

* have a working docker installation: https://www.docker.com/community-edition
* cd docker
* docker build . -t miker/jupyter
* cd ..
* command: `./startJupyter.sh`
* the startjupyter command returns a *docker_id*, copy it.
* command: `docker exec -it docker_id bash`

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

## Other

* AAVSO VSX catalog can be downloaded here: http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=B%2Fvsx%2Fversions%2F2018-02-26&target=brief&

## TODO

- new ensemble comparison star calculation?
- aperture is now fixed at 2, do some experiments on 1 star to measure effect on error bars and results
- check error column for any stars having error bars > 1%

## References

* https://books.google.ae/books?id=g_K3-bQ8lTUC&pg=PA236&lpg=PA236&dq=munipack+aperture&source=bl&ots=P4BKKI25HG&sig=Lj9Kg6EZi2pwKZXK5Hk5_B4qcVg&hl=en&sa=X&ved=0ahUKEwj4qr6o4evZAhWILcAKHdtnAqYQ6AEIRDAE#v=onepage&q&f=false
* https://www.aavso.org/how-report-new-variable-star-discoveries
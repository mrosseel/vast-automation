# munipack-automation

Docker is used to construct a light-weight virtual machine linux containing all necessary dependencies.
Once you're in this VM, all python commands can be run.

## Starting Docker

* have a working docker installation: https://www.docker.com/community-edition
* cd docker
* docker build . -t miker/jupyter
* cd ..
* command: `./startJupyter.sh`
* the startjupyter command returns a docker_id, copy it.
* command: `docker exec -it docker_id bash`

##Command line usage

###Init settings

* take one reference frame and calculate a fits header using http://Astrometry.net
* edit init.py to set all correct directories and values

###Run

* command: `python do_muniwin.py`

###File overview

* init.py : directory settings, stars to chart settings, ...
* do_muniwin.py : start all
* do_charts.py : only do charts
* do_upsilon.py : only do machine learning detection
* do_profile.py : do performance profiling on the app (not sure if working)

##Jupyter usage (deprecated)
* command: `jupyter notebook list`
* open browser with the provided url


## TODO

- calculate star position
- new ensemble comparison star calculation?
- capture upsilon extra data (minimum = period) to identify the main star?
- check error column for any stars having error bars > 1%

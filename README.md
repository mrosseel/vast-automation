# munipack-automation

## Starting jupyter

CMD: ./startJupyter.sh
get docker_id
CMD: docker exec -it docker_id bash
CMD: jupyter notebook list
open browser with the provided url

## Running python stuff

edit init.py to set all correct directories and values
python do_muniwin.py

* init.py : directory settings, stars to chart settings, ...
* do_muniwin.py : start all 
* do_charts.py : only do charts
* do_upsilon.py : only do machine learning detection
* do_profile.py : do performance profiling on the app (not sure if working)

TODO:
- capture upsilon extra data (minimum = period) to identify the main star
- check error column for any stars having error bars > 1%

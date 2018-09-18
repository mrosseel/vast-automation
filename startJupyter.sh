#!/bin/sh
docker stop variables
docker rm variables
docker run --privileged --name variables --user root -e GRANT_SUDO=yes -v $(pwd):/home/jovyan/work -d -p 8888:8888 miker/jupyter start.sh jupyter lab
# --NotebookApp.iopub_data_rate_limit=10000000 --NotebookApp.rate_limit_window=10.0
docker exec -it variables bash

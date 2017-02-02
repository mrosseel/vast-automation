docker run --user root -e GRANT_SUDO=yes -v $(pwd):/home/jovyan/work -d -p 8888:8888 miker/jupyter

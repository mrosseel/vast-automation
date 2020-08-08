FROM jupyter/scipy-notebook:7a0c7325e470
MAINTAINER Mike Rosseel version: 0.1
ENV CMUNIPACK_VERSION 2.1.27

USER root
RUN apt update && apt install -yq wcslib-dev libcfitsio-dev expat gtk2.0 gcc openssl python3-pip libfftw3-dev less nano  && \
# gstreamer-1.0 needed by cmunipack, some things could probably be removed
    apt install -yq libgstreamer1.0-0 gstreamer1.0-plugins-base gstreamer1.0-plugins-good gstreamer1.0-plugins-bad gstreamer1.0-plugins-ugly gstreamer1.0-libav gstreamer1.0-doc gstreamer1.0-tools gstreamer1.0-x gstreamer1.0-alsa gstreamer1.0-gl gstreamer1.0-gtk3 gstreamer1.0-qt5 gstreamer1.0-pulseaudio  && \
    apt install -yq build-essential gfortran g++ libx11-dev libxi-dev libxmu-dev libpng-dev curl wget bc iputils-ping && \
    wget https://github.com/kirxkirx/vast/archive/vast-1.0rc84.tar.gz && tar xvfz vast-1.0rc84.tar.gz && \
    mv vast-vast-1.0rc84 vast-1.0rc84
RUN pip install pipenv
WORKDIR /home/jovyan/vast-1.0rc84
RUN make
USER $NB_USER
WORKDIR /home/jovyan
RUN git clone https://github.com/mrosseel/vast-automation.git
RUN git clone https://github.com/dwkim78/upsilon.git  && \
    mkdir -p .jupyter
WORKDIR /home/jovyan/vast-automation
RUN pipenv install --system --deploy --ignore-pipfile
COPY jupyter_notebook_config.json .jupyter/jupyter_notebook_config.json
USER root
RUN chmod 777 .jupyter/jupyter_notebook_config.json
WORKDIR /home/jovyan/upsilon
RUN python3 setup.py install
WORKDIR /home/jovyan/work
EXPOSE 80
USER $NB_USER

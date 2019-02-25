
FROM jupyter/scipy-notebook:83ed2c63671f

MAINTAINER Mike Rosseel version: 0.1

ENV CMUNIPACK_VERSION 2.1.25

USER root
RUN apt-get update && apt-get install -yq wcslib-dev libcfitsio-dev expat gtk2.0 gcc openssl python3-pip libfftw3-dev less nano  && \
# gstreamer-1.0 needed by cmunipack, some things could probably be removed
    apt-get install -yq libgstreamer1.0-0 gstreamer1.0-plugins-base gstreamer1.0-plugins-good gstreamer1.0-plugins-bad gstreamer1.0-plugins-ugly gstreamer1.0-libav gstreamer1.0-doc gstreamer1.0-tools gstreamer1.0-x gstreamer1.0-alsa gstreamer1.0-gl gstreamer1.0-gtk3 gstreamer1.0-qt5 gstreamer1.0-pulseaudio  && \
    wget https://sourceforge.net/projects/c-munipack/files/C-Munipack%202.1%20Stable/2.1.25/cmunipack-2.1.25.tar.gz/download -O cmunipack.tar.gz  && \
    tar xvfz cmunipack.tar.gz
WORKDIR /home/jovyan/cmunipack-2.1.25
RUN chmod +x ./configure && ./configure && make && make install && ldconfig
USER $NB_USER
WORKDIR /home/jovyan
RUN pip install --upgrade pip && \
    pip install numpy && \
    pip install astropy pyreadline pyfftw tqdm photutils gatspy astroquery
# needs at least scikit 0.18.1
RUN conda update -y scikit-learn matplotlib && \
#RUN git clone https://github.com/mrosseel/upsilon.git
    git clone https://github.com/dwkim78/upsilon.git  && \
    mkdir -p .jupyter
COPY jupyter_notebook_config.json .jupyter/jupyter_notebook_config.json
USER root
RUN chmod 777 .jupyter/jupyter_notebook_config.json
WORKDIR /home/jovyan/upsilon
RUN python3 setup.py install
WORKDIR /home/jovyan/work
USER $NB_USER

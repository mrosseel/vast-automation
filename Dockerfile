
FROM jupyter/scipy-notebook

MAINTAINER Mike Rosseel version: 0.1

ENV CMUNIPACK_VERSION 2.1.15

USER root
RUN apt-get update && apt-get install -yq gstreamer1.0 wcslib-dev libcfitsio3-dev expat gtk2.0 gcc python3-pip libfftw3-dev
RUN wget https://sourceforge.net/projects/c-munipack/files/C-Munipack%202.1%20Stable/2.1.15/cmunipack-2.1.15.tar.gz/download -O cmunipack.tar.gz
RUN tar xvfz cmunipack.tar.gz
RUN cd cmunipack-2.1.15 && chmod +x ./configure && ./configure && make && make install && ldconfig
USER $NB_USER
#RUN pip2 install astropy pyfftw pyreadline &&
RUN pip3 install astropy pyreadline pyfftw tqdm
# needs at least scikit 0.18.1
RUN pip3 install --upgrade scikit-learn
RUN git clone https://github.com/mrosseel/upsilon.git
USER root
RUN cd upsilon && python3 setup.py install
USER $NB_USER

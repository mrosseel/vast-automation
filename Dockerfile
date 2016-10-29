
FROM jupyter/scipy-notebook

MAINTAINER Mike Rosseel version: 0.1

ENV CMUNIPACK_VERSION 2.1.14

USER root
RUN apt-get update && apt-get install -yq gstreamer1.0 wcslib-dev libcfitsio3-dev expat gtk2.0
RUN wget https://sourceforge.net/projects/c-munipack/files/C-Munipack%202.1%20Stable/2.1.14/cmunipack-2.1.14.tar.gz/download -O cmunipack-2.1.14.tar.gz
RUN tar xvfz cmunipack-2.1.14.tar.gz
RUN cd cmunipack-2.1.14 && chmod +x ./configure && ./configure && make && make install && ldconfig
USER $NB_USER
RUN pip2 install astropy upsilon pyreadline && pip3 install astropy upsilon pyreadline

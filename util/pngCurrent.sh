#!/bin/sh
cd current/converted_fits
fits2bitmap --stretch log *.fts
mv *.png ../animation
cd -

#!/bin/sh
wget --no-clobber --convert-links --random-wait -r -p --level 1 -E -e robots=off -U mozilla $1
#myurl=$1
#paths=${myurl#*/}
#lastpath=${}
# ln -s ./posts/gaia19bld_final/gaia19bld_final/index.html index.html

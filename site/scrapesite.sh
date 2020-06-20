#!/bin/sh
wget --no-clobber --convert-links --random-wait -r -p --level 1 -E -e robots=off -U mozilla $1

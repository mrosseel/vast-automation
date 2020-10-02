#!/bin/sh
# gets an exact copy of one hugo post. Afterwards, replace the top-level index.html with a symlink to the post index.html
wget --no-clobber --convert-links --random-wait -r -p --level 1 -E -e robots=off -U mozilla $1
python3 process_site.py $1
# echo "Now run find . -type f -name "*.png" -exec optipng -o5 {} \;"


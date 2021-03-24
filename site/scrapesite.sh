#!/bin/sh
# gets an exact copy of one hugo post. Afterwards, replace the top-level index.html with a symlink to the post index.html
# usage: ./scrapesite.sh http://ns313855.ip-188-165-231.eu/posts/wwcra2013-2017_ref1893/wwcra2013-2017_ref1893/
wget --no-clobber --convert-links -r -p --level 1 -E -e robots=off -U mozilla $1
python3 process_site.py $1
echo 'Now run find . -type f -name "*.png" -exec optipng -o5 {} \;'


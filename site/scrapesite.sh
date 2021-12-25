#!/bin/sh
# gets an exact copy of one hugo post. Afterwards, replace the top-level index.html with a symlink to the post index.html
#
# usage: ./scrapesite.sh http://ns313855.ip-188-165-231.eu/posts/wwcra2013-2017_ref1893/wwcra2013-2017_ref1893/
#
# After the png optimisation has run, you can rename the directories to avoid all the repetition. This needs surgery in index.html:
# - delete the star specific index.html, it's not needed
# - If you collapse the post directories to just one folder, rename ../../../images to ../../images in the post/.../index.html
#
wget --no-clobber --convert-links -r -p --level 1 -E -e robots=off -U mozilla $1
python3 process_site.py $1
echo 'Now run find . -type f -name "*.png" -exec optipng -o5 {} \;'


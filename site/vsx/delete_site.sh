#!/bin/sh
echo "Deleting posts and images for $@"
for dir in "$@"
do
    rm -Rf ./content/posts/$dir/
    rm -Rf ./static/images/$dir/
done

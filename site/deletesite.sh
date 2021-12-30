#/bin/sh
rm -Rf vsx/public/images/$1
rm -Rf vsx/public/posts/$1
rm -Rf vsx/content/posts/$1
vi vsx/public/index.xml

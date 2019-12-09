mv vsx.dat vsx.dat.old
curl --head --silent http://cdsarc.u-strasbg.fr/ftp/B/vsx/vsx.dat.gz | grep "Last-Modified" > vsx_last_modified.txt
curl -O http://cdsarc.u-strasbg.fr/ftp/B/vsx/vsx.dat.gz
gunzip vsx.dat.gz
python ./src/vsx_pickle.py vsx.dat

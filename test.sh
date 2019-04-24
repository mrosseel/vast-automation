#!/bin/bash
INTDIR=./inputfiles/integrationtest
find $INTDIR -maxdepth 1 -type f -exec rm -v {} \;
rm -Rf $INTDIR/results 
rm -Rf $INTDIR/conv*
rm -Rf $INTDIR/matched_photometry
rm -Rf $INTDIR/photometry
rm $INTDIR/*
cp $INTDIR/config_files/* $INTDIR
python src/set_reference_frame.py $INTDIR/fits/WWCrA#30V_000184532_FLAT.fit $INTDIR
timestamp=$(date +'%m-%d-%Y-%H_%M')
python -u ./src/do_muniwin.py -d $INTDIR --nowait --vsx --upsilon |& tee $INTDIR/test-$timestamp.log

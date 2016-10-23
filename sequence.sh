#!/bin/sh
cd inputfiles/WWCra
muniphot *.fit
munimatch --verbose frame001.pht frame*.pht
munilist -a 2 --object 1 -v 2 -c 3 out.txt match*.pht
munifind -a 2 -c 0 out2.txt match00*
cd ../..

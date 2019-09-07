#!/bin/bash
timestamp=$(date +'%m-%d-%Y-%H_%M')
python -u ./src/do_vast.py $* |& tee munilog-$timestamp.log

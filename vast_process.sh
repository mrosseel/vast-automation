#!/bin/bash
timestamp=$(date +'%m-%d-%Y-%H_%M')
python -u ./src/cli_vast.py $* 2>&1 | tee vastautomation-$timestamp.log

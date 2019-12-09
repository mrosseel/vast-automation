#!/bin/sh
ffmpeg -framerate 30 -i kout%06d.png -c:v libx264 -r 30 out.mp4
#!/usr/bin/env bash

ffmpeg -r 60 -f image2 -i ./build/anim_%05d.png -vf tpad=stop_mode=clone:stop_duration=4 -vcodec libx264 -crf 25 -pix_fmt yuv420p output.mp4

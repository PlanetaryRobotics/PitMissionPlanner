#!/usr/bin/env bash

ffmpeg -r 60 -f image2 -i $1/anim_%05d.png -vf tpad=stop_mode=clone:stop_duration=4 -vcodec libx264 -crf 25 -pix_fmt yuv420p $1/animation.mp4

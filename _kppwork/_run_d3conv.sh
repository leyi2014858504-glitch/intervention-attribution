#!/bin/bash
# 256-run n=8 game: 289-vs-288 daily-mean convention quantified end to end.
cd /mnt/d/atmos
python3 _kppwork/_d3_conv.py > _kppwork/_d3_conv.log 2>&1
echo "D3CONV EXIT $?"
tail -25 _kppwork/_d3_conv.log
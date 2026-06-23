#!/bin/bash
# 16^4 confirmation of the other-3 winner. Detached runner:
#   setsid bash search/confirm_launch.sh </dev/null >/dev/null 2>&1 &
cd ~/perfect_char || exit 1
while pgrep -x wilson_flow >/dev/null 2>&1; do sleep 15; done   # let any flows finish
rm -f runs/confirm.DONE
bash search/run_points.sh runs/confirm runs/confirm.pts 4 16 16 60 40 1500 8 4 0.04 9 > runs/confirm.log 2>&1
python3 search/analyze_scan.py runs/confirm >> runs/confirm.log 2>&1
touch runs/confirm.DONE

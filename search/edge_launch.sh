#!/bin/bash
# Edge-screen launcher for lenore. Waits for any in-flight Wilson flows (the dense
# scan) to finish, then runs the new-character edge screen + analysis. Meant to be
# started detached: `setsid bash search/edge_launch.sh </dev/null >/dev/null 2>&1 &`
cd ~/perfect_char || exit 1
# let the in-flight (dense) flows finish before we take the cores
while [ "$(pgrep -c wilson_flow 2>/dev/null || echo 0)" -gt 0 ]; do sleep 15; done
rm -f runs/edge.DONE
bash search/run_points.sh runs/edge runs/edge_screen.pts 4 12 12 60 40 1500 8 4 > runs/edge.log 2>&1
python3 search/analyze_scan.py runs/edge >> runs/edge.log 2>&1
touch runs/edge.DONE

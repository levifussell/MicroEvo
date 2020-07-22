#!/bin/sh

python auto_ep_slurm.py --save-location /lustre/answain/MicroEvo/results/ --kill-radius 0 10 --inhibit-radius 0 10 --grow-radius 0 10 --mutation 0 0.2 10 --kill-margin 0 0.5 10 --inhib-margin 0 0.5 10 --versions 4 5 --num-jobs 128

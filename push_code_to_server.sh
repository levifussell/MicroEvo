#!/bin/sh

rsync -azv --progress --filter="- /data" --filter="- *.csv" --filter="- .*" --filter="- *.ipynb" . answain@login.deepthought2.umd.edu:MicroEvo/

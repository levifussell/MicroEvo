#!/bin/sh

file_location="data/PLOT3/data/"
kill_radius=(0 3)
inhib_radius=(0 4)
grow_radius=(0 5)
mutation=(0.0 0.1 2)
kill_margin=(0.0 0.3 2)
inhib_margin=(0.0 0.3 2)
version=3

./a.out $file_location ${kill_radius[0]} ${kill_radius[1]} ${inhib_radius[0]} ${inhib_radius[1]} ${grow_radius[0]} ${grow_radius[1]} ${mutation[0]} ${mutation[1]} ${mutation[2]} ${kill_margin[0]} ${kill_margin[1]} ${kill_margin[2]} ${inhib_margin[0]} ${inhib_margin[1]} ${inhib_margin[2]} ${version}



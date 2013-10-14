#!/bin/bash
gnuplot -e "unset colorbox; set xrange[0:62]; set yrange[0:62]; set pm3d map; set pm3d interpolate 0,0; set terminal wxt size 1000,700; splot '$1' matrix;"
#gnuplot -e "unset colorbox; set pm3d map; set pm3d interpolate 0,0; set terminal wxt size 1000,700; set cntrparam levels 20; set contour; splot '$1' matrix;"

#!/usr/bin/gnuplot

set xrange[0:0.01]
set yrange[0:0.01]

set zrange[0:2000]
set cbrange[0:2000]

do for[i=1:100]{splot 'fichier/T'.i.'.dat' u 1:2:3 w pm3d title 'T'.i.'.dat';pause 0.1}

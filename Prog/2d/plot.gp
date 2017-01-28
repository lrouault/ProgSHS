#!/usr/bin/gnuplot

set xrange[0:0.01]
set yrange[0:0.01]

set zrange[0:2000]
set cbrange[0:2000]

#do for[i=1:100]{
#splot 'fichier/T'.i.'P0.dat' u 1:2:3 w p title 'T'.i.'P0.dat',\
#'fichier/T'.i.'P1.dat' u 1:2:3 w p title 'T'.i.'P1.dat',\
#'fichier/T'.i.'P2.dat' u 1:2:3 w p title 'T'.i.'P2.dat',\
#'fichier/T'.i.'P3.dat' u 1:2:3 w p title 'T'.i.'P3.dat';
#pause 0.1}

set pm3d map
set xlabel 'Nx'
set ylabel 'Ny'
splot 'fichier/T0P0.dat' u 1:2:3 w l title 'P0',\
'fichier/T0P1.dat' u 1:2:3 w l title 'P1',\
'fichier/T0P2.dat' u 1:2:3 w l title 'P2',\
'fichier/T0P3.dat' u 1:2:3 w l title 'P3',\
'fichier/T0P4.dat' u 1:2:3 w l title 'P4',\
'fichier/T0P5.dat' u 1:2:3 w l title 'P5'

pause -1

set xrange[0:0.01]
set yrange[0:0.01]

set zrange[0:1]
set cbrange[0:1]


i=a
splot 'T'.i.'.dat' u 1:2:4 w pm3d title 'T'.i.'.dat'
pause 0.1
a=a+1
if (i<b) reread

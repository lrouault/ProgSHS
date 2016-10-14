set xrange[0:0.5]
set yrange[0:0.1]
set zrange[0:3000]


i=a
splot 'fichier/T'.i.'.dat' w pm3d title 'T'.i.'.dat'
pause 0.2
a=a+1
if (i<b) reread

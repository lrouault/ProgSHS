set xrange[0:202]
set yrange[0:202]
set zrange[0:1900]
set pm3d


do for [i=0:100]{splot 'fichier/T'.i.'.dat', pause=2}

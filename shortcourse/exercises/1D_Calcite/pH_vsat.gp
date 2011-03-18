set xtics font "Helvetica, 20"
set x2tics font "Helvetica, 20"
set ytics font "Helvetica, 20"
set xrange[4.8:8.7]
set x2range[0:1.2]
set yrange[0:10]
set xlabel 'pH' font "Helvetica,24"
set x2label 'Saturation [-]' font "Helvetica,24"
set ylabel 'Z [m]' font "Helvetica,24"
set key bottom right title 'Years' font "Helvetica,24"
plot \
'pflotran-001.tec' u 6:3 t '5' w l lw 2 lt 1,\
'pflotran-002.tec' u 6:3 t '10' w l lw 2 lt 2,\
'pflotran-003.tec' u 6:3 t '15' w l lw 2 lt 4,\
'pflotran-004.tec' u 6:3 t '20' w l lw 2 lt 5,\
'pflotran-005.tec' u 6:3 t '25' w l lw 2 lt 6,\
'pflotran-005.tec' u 5:3 t 'sat' axes x2y1 w l lw 3 lt 3


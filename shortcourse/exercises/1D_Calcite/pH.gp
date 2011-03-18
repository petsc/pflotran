set xtics font "Helvetica, 20"
set ytics font "Helvetica, 20"
set yrange[4.8:8.2]
set xlabel 'X [m]' font "Helvetica,24"
set ylabel 'pH' font "Helvetica,24"
set key bottom right title 'Years' font "Helvetica,24"
plot \
'pflotran-001.tec' u 1:4 t 'OS: 5' w l lw 2 lt 1,\
'pflotran-002.tec' u 1:4 t 'OS:10' w l lw 2 lt 2,\
'pflotran-003.tec' u 1:4 t 'OS:15' w l lw 2 lt 4,\
'pflotran-004.tec' u 1:4 t 'OS:20' w l lw 2 lt 5,\
'pflotran-005.tec' u 1:4 t 'OS:25' w l lw 2 lt 6


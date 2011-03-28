set xtics font "Helvetica, 20"
set ytics font "Helvetica, 20"
set y2tics font "Helvetica, 20"
set yrange[1.e-5:.1]
set y2range[0:1.2e-5]
set log y
set xlabel 'X [m]' font "Helvetica,24"
set ylabel 'Concentration [M]' font "Helvetica,24"
set y2label 'Volume Fraction [-]' font "Helvetica,24"
set key top left title '25 Years' font "Helvetica,24"
plot \
'pflotran-005.tec' u 1:5 t 'H+' w l lw 2 lt 1,\
'pflotran-005.tec' u 1:6 t 'HCO3-' w l lw 2 lt 2,\
'pflotran-005.tec' u 1:7 t 'Ca++' w l lw 2 lt 4,\
'pflotran-005.tec' u 1:8 t 'Calcite' axes x1y2 w l lw 2 lt 5


reset
set font "Helvetica, 20"
set lmargin 22
set bmargin 5
set tmargin 2
set rmargin 24
set xtics font "Helvetica, 20" offset character 0,-.5
set ytics font "Helvetica, 20"
set y2tics font "Helvetica, 20"
set format y2 "%.1g"
set yrange[1.e-5:.1]
set y2range[0:1.2e-5]
set log y
set xlabel 'X [m]' font "Helvetica,24" offset character 0,-1.25 
set ylabel 'Concentration [M]' font "Helvetica,24" offset character -8,0
set key inside top left Left title '25 Years' font "Helvetica,18" spacing 1.5
plot \
'pflotran-005.tec' u 1:5 t 'H+' w l lw 2 lt 1,\
'pflotran-005.tec' u 1:6 t 'HCO3-' w l lw 2 lt 2,\
'pflotran-005.tec' u 1:7 t 'Ca++' w l lw 2 lt 4,\
'pflotran-005.tec' u 1:8 t 'Calcite' axes x1y2 w l lw 2 lt 5


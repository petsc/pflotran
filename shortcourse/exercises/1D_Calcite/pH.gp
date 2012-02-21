reset
set font "Helvetica, 20"
set lmargin 18
set bmargin 5
set tmargin 2
set xtics font "Helvetica, 20"
set ytics font "Helvetica, 20"
set yrange[4.8:8.2]
set xlabel 'X [m]' font "Helvetica,24" offset character 0,-1.25
set ylabel 'pH' font "Helvetica,24" offset character -8,0
set key inside bottom right Left title 'Years' font "Helvetica,18" spacing 1.5
plot \
'pflotran-001.tec' u 1:4 t '5' w l lw 2 lt 1,\
'pflotran-002.tec' u 1:4 t '10' w l lw 2 lt 2,\
'pflotran-003.tec' u 1:4 t '15' w l lw 2 lt 4,\
'pflotran-004.tec' u 1:4 t '20' w l lw 2 lt 5,\
'pflotran-005.tec' u 1:4 t '25' w l lw 2 lt 6


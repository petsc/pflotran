reset
set size 0.95, 0.95
set bmargin 5
set lmargin 16
unset key
set xtics font "Helvetica, 20"
set ytics font "Helvetica, 20"
set xrange[0:3.5]
set yrange[0:0.3]
set ylabel 'Saturation [-]' offset -5,0 font "Helvetica,24"
set xlabel 'Time [y]' offset 0,-1 font "Helvetica,24" 
plot \
'observation-0.tec' u 1:3 t 'sat' axes x1y1 w l lw 3 lt 3


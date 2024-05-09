set term pngcairo size 800,800 enhanced
set output 're-atan.png'
set size ratio 0.8
set xrange [0:15]
set yrange [0:15]

w = 2

set title "Cylinder" font ",16"
set xlabel "time" font ",14"
set ylabel "CD" font ",14"
set tics font ",12"


plot 'log' u 2:($3==0?$8:NaN) w l t 'Re = 1' lw 3, 'log' u 2:($3==1?$8:NaN) w l t 'Re = 2' lw 3, 'log' u 2:($3==2?$8:NaN) w l t 'Re = 5' lw 3, 'log' u 2:($3==3?$8:NaN) w l t 'Re = 10' lw 3, 'log' u 2:($3==4?$8:NaN) w l t 'Re = 20' lw 3, 'log' u 2:($3==5?$8:NaN) w l t 'Re = 40' lw 3

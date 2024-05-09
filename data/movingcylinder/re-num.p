set term pngcairo size 800,800 enhanced
set output 're-2000.png'
set size ratio 0.8
set xrange [0:15]
set yrange [0:275]

set title "Re = 2000" font ",16"
set xlabel "time" font ",14"
set ylabel "Enstrophy" font ",14"
set tics font ",12"


plot 'loga4' u 2:($3==0?$10:NaN) w l t 'IBM smeared' lw 3, 'loga3' u 2:($3==0?$10:NaN) w l t 'IBM' lw 3, 'loga5' u 2:($3==0?$10:NaN) w l t 'very simple IBM' lw 3, 'loga1' u 2:($3==0?$10:NaN) w l t 'embed' lw 3


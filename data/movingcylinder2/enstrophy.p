set term pngcairo size 800,800 enhanced
set output 'enstrophy-re-40.png'
set size ratio 0.8
set xrange [0:15]
set yrange [0:15]

set title "Cylinder" font ",16"
set xlabel "time" font ",14"
set ylabel "CD" font ",14"
set tics font ",12"

a = 2      # time
b = '$8'  # 8 = CD, 10 = enstrophy
j = 0      # j (0: Re = 1, 1: Re = 2, etc.)


plot 'log' u a:($3==j?$8:NaN) w l t 'incorrect smoothening' lw 3, \
'log1' u a:($3==j?$8:NaN) w l t 'correct smoothening' lw 3, \
'log2' u a:($3==j?$8:NaN) w l t 'embed' lw 3
#'log3' u a:($3==j?$8:NaN) w l t 'no smoothening' lw 3




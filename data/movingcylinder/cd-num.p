set term pngcairo size 900,900 enhanced
set output 'cd-2000.png'
set size ratio 0.8
set xrange [0:15]
set yrange [0:1]

set title "Re = 2000" font ",16"
set xlabel "time" font ",14"
set ylabel "CD" font ",14"
set tics font ",12"


plot 'loga4' u 2:($3==0?$9:NaN) t 'IBM smeared', 'loga3' u 2:($3==0?$9:NaN) t 'IBM' , 'loga5' u 2:($3==0?$9:NaN) t 'very simple IBM' , 'loga1' u 2:($3==0?$9:NaN) t 'embed'


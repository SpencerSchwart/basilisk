set term pngcairo size 800,800 enhanced
set output 'vlogcd.png'
set size ratio 0.8
set xrange [0:20]
set yrange [0:3]

set title "CD vs Time" font ",16"
set xlabel "Time" font ",14"
set ylabel "CD" font ",14"
set tics font ",12"
set pointsize 1.5

plot 'vlog' u 2:11 t 'SPM' with points pointtype 7 lc "green", \
   'vlog4' u 2:11 t 'no SPM' with points pointtype 1 lc "blue", \
    'vlog3'u 2:8 t 'embed' with points pointtype 13 lc "red", \
      'vlogtest' u 2:12 t 'full SPM' w l lc "orange"



set term pngcairo size 800,800 enhanced
set output 'vprof.png'
set size ratio 0.8
set xrange [2.5:3.5]
set yrange [-0.1:0.2]

set title "Centerline Velocity" font ",16"
set xlabel "y" font ",14"
set ylabel "U.x" font ",14"
set tics font ",12"
set pointsize 1.5
set key center top

plot 'vprofy' u 3:4 t 'SPM' with points pointtype 7 lc "green", \
   'vprofy4' u 3:4 t 'no SPM' with points pointtype 1 lc "blue", \
    'vprofy3'u 3:4 t 'embed' with points pointtype 13 lc "red", \
      'vprofytest' u 3:4 t 'full SPM' w p pointtype 24 lc "orange"



set term pngcairo size 800,800 enhanced
set output 'cd-v-re.png'
set size ratio 0.8
set xrange [0:50]
set yrange [0:15]

set title "CD vs Re" font ",16"
set xlabel "Re" font ",14"
set ylabel "CD" font ",14"
set tics font ",12"
set pointsize 1.5

a= 2      # 2 = Re
# 3 = CD, 4 = enstrophy

# Enstrophy
#plot 'log-ss' u a:($1==1?$4:NaN) t 'smoothening IBM' with points pointtype 7 lc "green", \
#   'log-ss' u a:($1==3?$4:NaN) t 'no smoothening IBM' with points pointtype 1 lc "blue", \
#      'log-ss' u a:($1==2?$4:NaN) t 'embed' with points pointtype 9 lc "orange", \
 #     'log-ss' u a:($1==4?$4:NaN) t 'improved IBM' with points pointtype 13 lc "red"
#        'log-ss' u a:($1==5?$4:NaN) t 'very simple IBM' with point pointtype 24 lc "cyan"

# CD
plot 'default-data' u 1:2 t 'Paper' with points pointtype 5 lc "purple", \
'log-ss' u a:($1==1?$3:NaN) t 'smoothening IBM' with points pointtype 7 lc "green", \
   'log-ss' u a:($1==3?$3:NaN) t 'no smoothening IBM' with points pointtype 1 lc "blue", \
      'log-ss' u a:($1==2?$3:NaN) t 'embed' with points pointtype 9 lc "orange", \
      'log-ss' u a:($1==4?$3:NaN) t 'improved IBM' with points pointtype 13 lc "red", \
#        'log-ss' u a:($1==5?$3:NaN) t 'very simple IBM' with point pointtype 24 lc "cyan"







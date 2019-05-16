#set terminal windows
set terminal png font arial 14 size 1600,600
set output 'T(t).png'
#----------------------------------------------------
set title "Зависимость температуры от времени"
set grid
set key top right

# X
set xlabel "t, мсек" font ",18"
#set xrange [9:19] 

# Y
set ylabel "T, K" font ",18"
set format y "%2.2f"
#set yrange [0:3000] 
#set format y '%.0s✕10^{%S}' 


set pointsize 2
set style line 5 lt rgb "cyan" lw 3 pt 6
set style function linespoints
set style line 1 linetype 1 linecolor rgb "black"


plot "../postProcessing/P(t)/0/T" using ($1*1000.0):($2) title "" with lines linetype 1 linecolor "red" linewidth 2

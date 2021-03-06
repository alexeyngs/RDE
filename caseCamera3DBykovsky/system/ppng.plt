set terminal png size 1600,600 font "arial,18"
set output 'p(t).png'
#----------------------------------------------------
set title "Зависимость давления от времени"
set grid
set key top right

# X
set xlabel "t, мсек" font ",18"
#set xrange [9:19] 

# Y
set ylabel "P, атм" font ",18"
set format y "%2.2f"
#set yrange [0:3000] 
#set format y '%.0s✕10^{%S}' 


set pointsize 2
set style line 5 lt rgb "cyan" lw 3 pt 6
set style function linespoints
set style line 1 linetype 1 linecolor rgb "black"


plot "../postProcessing/P(t)/0/p" using ($1*1000.0):($2/101325.0) title "" with lines linetype 1 linecolor "black" linewidth 1

#pause -1
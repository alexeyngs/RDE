set terminal png font arial 14 size 800,600
set output 'P(t).png'
#----------------------------------------------------
set title "Зависимость давления от времени"
set grid

set xlabel "t, сек" font	",18"	
set ylabel "P" font	",18"	

set pointsize 2
set style line 10 linetype 1 linecolor rgb "black"
set style function linespoints


plot "p" with lines lt 1 linecolor "black"

#
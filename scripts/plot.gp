reset
set xlabel 'F(n)'
set ylabel 'time (ns)'
set title 'Fibonacci execution time'
set term png enhanced font 'Verdana,10'
set output 'plot_output.png'
plot [2:100][:] \
'plot_input' using 1:2 with linespoints linewidth 2 title "recursion w/ cache",\
'' using 1:3 with linespoints linewidth 2 title "fast doubling w/ clz",\
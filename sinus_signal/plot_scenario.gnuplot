set datafile separator ',' # use comma as separator in data file, it is default separator for CSV files
set key outside # legend outside of plot
set key autotitle columnhead # use the first line as title
set xlabel 'Time' # label for the X axis
set ylabel "Data" # label for the Y axis
set grid # use grid with plot
plot for [i=2:*] "data.csv" using 1:i with linespoints pt 7 # plot all columns from "data.csv" file starting from second column, use first column as X axis data

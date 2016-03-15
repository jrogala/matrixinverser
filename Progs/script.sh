#! /bin/sh
g++ -Wall benchmark.cc matrix.cc -o benchmark
./benchmark
gnuplot -e "plot 'data_gauss.dat' using 1:2 title 'Gauss time'"

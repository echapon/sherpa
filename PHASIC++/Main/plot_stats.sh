#!/bin/bash

test -z "$1" && exit 1;

rdb=$1;
n=0;
for i in $(unzip -l $rdb '*Statistics.dat' | grep Statistics | awk '{print $NF}'); do
  (( ++n ));
  mkdir -p $(dirname $i);
  unzip -p $rdb $i > $i;
  if test -z "$plotcmd"; then plotcmd="plot ";
  else plotcmd=$plotcmd", "; fi;
  t=$(echo $i | sed -e 's|.*MC_._.__||g;s|/Statistics.dat||g;s|__QCD.*||g');
  fitcmd=$fitcmd"f$n(x) = a$n; fit f$n(x) '$i' using 5:(\$2):(2*\$3) via a$n;"
  plotcmd=$plotcmd"'$i' u 5:(abs(\$2)):(2*\$3) w yerr t '$t' lc $n lt $n"
  plotcmd=$plotcmd", abs(a$n) t sprintf(\"%3.4g pb\",a$n) lc $n lt $n";
done;

gnuplot <<EOF
set term postscript color;
set output 'stats_plot.ps';
set logscale xy;
set xlabel "number of points"
set ylabel "cross section [pb]"
set size 1.2, 0.9
set key outside
set key reverse
set key Left
$fitcmd
$plotcmd
EOF

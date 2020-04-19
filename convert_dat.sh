SCALE=$3
XMAX=`echo 3.22*${SCALE}|bc -l`
YMAX=`echo 1.2*${SCALE}|bc -l`

cat <<EOSCRIPT | gnuplot -
set term png size 1200,1200
set output "$2"
set style data dots
set xrange [0:${XMAX}]
set yrange [0:${YMAX}]
plot "$1" using 2:3
EOSCRIPT

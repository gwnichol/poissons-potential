# to use this from within gnuplot type: load 'graph.txt', but first change the directorty (click on the ChDir button and select yours)
reset

set terminal pngcairo enhanced font "Garamond,12" fontscale 1.0 size 1000,800 
set output ARG2 # we assume here that you used the ChDir button in wgnuplot to select your working directory...

set style line 100 lt 1 lc 0 lw 0.2 #redefine grid style
set grid ls 100

#set xtics 2 offset 0,-0.5
#set ytics 2 offset 1,0
#set ztics 2
#set mxtics 2
#set mytics 2
#set mztics 2
#set ticslevel 0.3
#set border 895

set xlabel 'Position x (cm)' offset 0,-1.5
set ylabel 'Position y (cm)' offset 3.5, 0
set zlabel 'Potential (V)' rotate

set style data pm3d
set palette defined (0 "#663c5c", 1 "#ff96e8", 2 "black", 3 "#96a0ff", 4 "#3c4066")

set view equal xy
set view ARG3,ARG4,1,1
set contour base
set cntrparam levels incremental  -10.0,0.25,10.0

unset key
splot ARG1 u 1:2:3 t "" 

reset
set term wxt
set out

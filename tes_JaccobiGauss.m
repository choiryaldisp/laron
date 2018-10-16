clear
close
clc

%matrix input
A=[1 2 -2; 1 1 1; 2 2 1];
bi=[7 2 5];
x0=[0,0,0];

%stop condition
epsilon=10^-5;
iter_max=25;

%Calculate Progress 
[ solG, all_solG, N_iterG ] = GaussSiedel( A, bi, x0, iter_max, epsilon );
[ solJ, all_solJ, N_iterJ ] = Jaccobi( A, bi, x0, iter_max, epsilon );
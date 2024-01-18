% demo_sugr_Point_2D: demo routine for generating 2D point
%
% Point_2D = structure 
%
% *           .h = spherically normalized homogeneous coordinates
% *           .Crr = half vector of reduced covariance matrix
% *           .type = 'Point_2D'
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de

%% initialize sugr
addpath(genpath('../../General-Functions'))
close all; clearvars;clc
sugr_INIT
plot_init;

%% Sorted according to number of arguments

%% (1) Euclidean point vector without CovM (==0)
in = [1,2]';
x = sugr_Point_2D(in);
sugr_show_Point_2D(x,'x');


%% (1) Homogeneous point vector without CovM (==0)
in = [1,2,3]';
x = sugr_Point_2D(in);
sugr_show_Point_2D(x,'x');

%% (1) Homogeneous vanishing point vector without CovM (==0)
in = [1,2,0]';
x = sugr_Point_2D(in);
sugr_show_Point_2D(x,'x');

%% (2) Euclidean point coordinates without CovM (==0)
in1 = 1;
in2 = 2;
x = sugr_Point_2D(in1,in2);
sugr_show_Point_2D(x,'x');

%% (2) Homogeneous point vector 
in1 = [1,2,2]';
in2 = [1,2,1;2,5,0;1,0,7]*0.01;
x = sugr_Point_2D(in1,in2);
sugr_show_Point_2D(x,'x');

%% (2) Homogeneous vanishing point vector
in1 = [3,4,0]';
in2 = [2,2,0;2,6,0;0,0,5];
x = sugr_Point_2D(in1,in2);
sugr_show_Point_2D(x,'x');

%% (3) Homogeneous point coordinatens without CM (==0)
in1 = 1;
in2 = 2;
in3 = 2;
x = sugr_Point_2D(in1,in2,in3);
sugr_show_Point_2D(x,'x');

%% (3) Homogeneous vanishing point coordinates without CM (==0)
in1 = 1;
in2 = 2;
in3 = 0;
x = sugr_Point_2D(in1,in2,in3);
sugr_show_Point_2D(x,'x');

%% (4) Euclidean point coordinates
in1 = 1;
in2 = 2;
in3 = 0.1;
in4 = 0.2;
x = sugr_Point_2D(in1,in2,in3,in4);
sugr_show_Point_2D(x,'x');

%% (4) Homogeneous point coordinates
in1 = 1;
in2 = 2;
in3 = 2;
in4 = [1,2,1;...
       2,5,0;...
       1,0,7]*0.0001;
x1 = sugr_Point_2D(in1,in2,in3,in4);
sugr_show_Point_2D(x,'x_1');
figure('Color','w')
sugr_plot_Point_2D(x1,'-k','-b',2,1);
title('Homogeneous uncertain point')
axis equal

%% (4) Homogeneous vanishing point coordinates
in1 = 1;
in2 = 2;
in3 = 0;
in4 = [2,2,1;2,5,0;1,0,7]*0.01;
x = sugr_Point_2D(in1,in2,in3,in4);
sugr_show_Point_2D(x,'x');


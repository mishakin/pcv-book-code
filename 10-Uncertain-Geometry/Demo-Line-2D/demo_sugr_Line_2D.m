% demo_sugr_Line_2D: demo routine for creating a 2D lines
%
% Line_2D = structure
%
% *           .h = spherically normalized homogeneous coordinates
% *           .Crr = reduced covariance matrix of l.h
% *           .type = 2
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de


% clear all
clearvars
close all
clc

addpath(genpath('../../General-Functions'))
%% initialize sugr
sugr_INIT

ss = plot_init;

%% (1) line through origin with Hessian parameter vector, CovM == 0
in = [atan2(4,3),0]';
l = sugr_Line_2D(in);
sugr_show_Line_2D(l,'l');

%% (1) line with Hessian parameters, CovM == 0
in = [3/5,4/5,-5]';
l = sugr_Line_2D(in);
sugr_show_Line_2D(l,'l');

%% (2) line with Hessian parameters,, CovM == 0
in1 = atan2(4,3);
in2 = 5;
l = sugr_Line_2D(in1,in2);
sugr_show_Line_2D(l,'l');

%% (2) line with Hessian parameter vector
in1 = [atan2(4,3),5]';
in2 = [0.0002 0.001;0.001 0.01]*0.01';
l = sugr_Line_2D(in1,in2);
sugr_show_Line_2D(l,'l');

figure('Color','w','Position',[50,100,0.45*ss(1),0.45*ss(1)]);
sugr_plot_Line_2D(l,'-b','-b',2,20);
title('line with Hessian parameter vector')
axis equal
% xlim([-4,10]);ylim([-1,10]);

%% (2) line with homogeneous coordinate vector
in1 = [0.1177,0.1569,-0.9806]';
in2 = 10^(-4)*[0.1714    0.0521    0.0289;...
    0.0521    0.0477    0.0139;...
    0.0289    0.0139    0.0057]';
l = sugr_Line_2D(in1,in2);
sugr_show_Line_2D(l,'l');

%% (2) line at infinity with homogeneous coordinate vector
in1 = [0 0 2]';
in2 = 10^(-4)*[0.1714    0.0521    0.0289;...
    0.0521    0.0477    0.0139;...
    0.0289    0.0139    0.0057]';
l = sugr_Line_2D(in1,in2);
sugr_show_Line_2D(l,'l');

%% (3) line with homogeneous coordinates, CovM == 0
in1 = 0.1177;
in2 = 0.1569;
in3 = -0.9806;
l = sugr_Line_2D(in1,in2,in3);
sugr_show_Line_2D(l,'l');

%% (3) line with Hessian parameters
in1 = atan2(4,3);
in2 = 5;
in3 = [0.0002 0.001;0.001 0.01]';
l = sugr_Line_2D(in1,in2,in3);
sugr_show_Line_2D(l,'l');

%% (4) line with homogeneous coordinates
in1 = 0.1177;
in2 = 0.1569;
in3 = -0.9806;
in4 = 10^(-4)*[0.1714    0.0521    0.0289;...
   0.0521    0.0477    0.0139;...
   0.0289    0.0139    0.0057]';
l = sugr_Line_2D(in1,in2,in3,in4);
sugr_show_Line_2D(l,'l');

%% (5) line with centroid coordinates, CpvM == 0
in1 = -1;
in2 = 7;
in3 = atan2(4,3);
in4 = 0.0141;
in5 = 0.0707;
l = sugr_Line_2D(in1,in2,in3,in4,in5);
sugr_show_Line_2D(l,'l');

figure('Color','w','Position',[0.55*ss(1),100,0.45*ss(1),0.45*ss(1)]);
sugr_plot_Line_2D(l,':r','-r',2,5);
title('line with centroid coordinates, CpvM == 0')
axis equal
% xlim([-4,8]);ylim([-1,8]);
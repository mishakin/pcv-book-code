% demo_construct_Point_2D
% demo for constructing points with SUGR
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de

clearvars
close all
clc

addpath(genpath('../../General-Functions'))
sugr_INIT
ss = plot_init;

%% intersection
l = sugr_Line_2D(0,1,3*pi/2,0.01,0.02);
sugr_show_Line_2D(l,'l');
m = sugr_Line_2D(1,0,pi/4,.01,0.01);
sugr_show_Line_2D(m,'m');
n = sugr_Line_2D([0,0,1]');

x = sugr_construct_intersection_Point_2D(l,m);
sugr_show_Point_2D(x,'x');
winkel_lm = sugr_angle_Line_2D(m,l);
% winkel_lm.a*180/pi
winkel_nm = sugr_angle_Line_2D(m,n);
fprintf('\nangle l-m = %5.3f°\tangle n-m = %5.3f°\n',winkel_lm.a*180/pi, winkel_nm.a*180/pi)

figure('name','intersection of two lines l and m','Color','w','Position',[150,0.28*ss(1),0.3*ss(1) 0.4*ss(2)])
hold on
sugr_plot_Line_2D(l,':r','-r',2,5);
sugr_plot_Line_2D(m,':b','-b',2,5);
sugr_plot_Point_2D(x,'.k','-k',2,5);
axis equal
title('intersection of two lines l and m')

% return

%% footpoint x
l = sugr_Line_2D(0,1,2*pi/3,0.01,0.01);
sugr_show_Line_2D(l,'l');
x = sugr_Point_2D(-1,1,0.01,0.01);
sugr_show_Point_2D(x,'x');

z = sugr_construct_footpoint_xl_Point_2D(x,l);

figure('name','footpoint x on line l','Color','w','Position',[0.4*ss(1),0.28*ss(1),0.3*ss(1) 0.4*ss(2)])
sugr_plot_Line_2D(l,'-r','-r',2,5);
sugr_plot_Point_2D(x,'-b','-b',2,5);
sugr_plot_Point_2D(z,'.k','-k',2,5);
axis equal
title('footpoint x on line l')
% 

%% midpoint of x and y
x = sugr_Point_2D(0.1*[-1;0],0.00001*[0.7, 1;1,2]);
sugr_show_Point_2D(x,'x');
y = sugr_Point_2D(0.1*[+1;0],0.00001*[0.7,-1;-1,2]);
sugr_show_Point_2D(y,'y');

z = sugr_construct_midpoint_Point_2D(x,y);
figure('name','midpoint of x and y','Color','w','Position',[0.3*ss(1),50,0.3*ss(1) 0.4*ss(2)])
sugr_plot_Point_2D(x,'.r','-r',2,5);
sugr_plot_Point_2D(y,'.b','-b',2,5);
sugr_plot_Point_2D(z,'.k','-y',2,5);
axis equal
title('midpoint and meanpoint of x and y')
% 

% mean point
zm = sugr_construct_mean_Point_2D(x,y);
sugr_show_Point_2D(zm,'zm')
sugr_plot_Point_2D(zm,'.k','-k',2,5);


%% footpoint x from O to l
l = sugr_Line_2D(0,1,2*pi/3,0.03,0.03);
sugr_show_Line_2D(l,'l');
x = sugr_Point_2D(0,0,10^(-3),10^(-3));
sugr_show_Point_2D(x,'x');

z = sugr_construct_footpoint_Ol_Point_2D(l);
figure('name','footpoint x from O to l','Color','w','Position',[0.65*ss(1),50,0.3*ss(1) 0.4*ss(2)])
sugr_plot_Line_2D(l,'-r','-r',2,5);
sugr_plot_Point_2D(x,'-b','-b',4,5);
sugr_plot_Point_2D(z,'.k','-k',2,5);
title('footpoint x from O to l')
axis equal
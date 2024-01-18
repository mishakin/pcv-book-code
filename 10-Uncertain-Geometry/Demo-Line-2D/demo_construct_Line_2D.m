%% test_construct_Line_2D: test routine for constructing 2D lines
%
% see PCV Sect. 10, Fig. 10.1
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de


close all
addpath(genpath('../../General-Functions'))

%% initialize sugr
sugr_INIT
ss = plot_init;

%% join and midline
% define 3 points
x = sugr_Point_2D([1,2]',0.00001*[2,2;2,6]);
y = sugr_Point_2D([-1,0]',0.00001*[2,-1;-1,2]);
u = sugr_Point_2D([0,2]',0.00001*[2,1;1,1]);

%% midpoint
z = sugr_construct_midpoint_Point_2D(x,y);

%%
% joining line
l = sugr_construct_join_Line_2D(x,y);
sugr_show_Line_2D(l,'l');
%%
% midline : check on commutativity 
m1 = sugr_construct_midline_Line_2D(x,y);
m2 = sugr_construct_midline_Line_2D(y,x);
%%
% parallel to origin
lp = sugr_construct_parallel_lx_Line_2D(l,u);

figure('Color','w','Position',[20 ss(2)/3 ss(1)/3 ss(2)/2])
hold on
xlim([-2,2]);
ylim([-1,3]);
sugr_plot_Point_2D(x,'.k','-b',2,40);
sugr_plot_Point_2D(y,'.k','-b',2,40);
sugr_plot_Point_2D(z,'.k','-b',2,40);
sugr_plot_Point_2D(u,'.k','-b',2,40);
sugr_plot_Line_2D(l,'-k','-r',2,40);
%sugr_plot_Line_2D(m1,'-k','-r',2,40);
%sugr_plot_Line_2D(lp,'-k','-r',2,40);
axis equal
title('Joining line ')

figure('Color','w','Position',[100+ss(1)/3  ss(2)/3 ss(1)/3 ss(2)/2])
hold on
xlim([-2,2]);
ylim([-1,3]);
sugr_plot_Point_2D(x,'.k','-b',2,40);
sugr_plot_Point_2D(y,'.k','-b',2,40);
sugr_plot_Point_2D(z,'.k','-b',2,40);
sugr_plot_Line_2D(l,'-k','-r',2,40);
sugr_plot_Line_2D(m2,'-k','-r',2,40);
axis equal
title('Midline ')

%% similar to Fig. 10.1
% define 2D points
x = sugr_Point_2D([-2,1]',   0.000001*[2,2;2,6]);
y = sugr_Point_2D([-1,1.5]', 0.000001*[2,-1;-1,2]);
z = sugr_Point_2D([1,1.3]',  0.000001*[2,-2;-2,6]);
u = sugr_Point_2D([1,-0.2]',0.000001*[2,1;1,2]);

l = sugr_construct_join_Line_2D(x,y);
m = sugr_construct_join_Line_2D(z,u);

v = sugr_construct_intersection_Point_2D(l,m);

figure('Color','w','Position',[50+ss(1)/6  50 ss(1)/3 ss(2)/2])
hold on
xlim([-2,2]);
ylim([-1,3]);
sugr_plot_Point_2D(x,'.k','-b',2,40);
sugr_plot_Point_2D(y,'.k','-b',2,40);
sugr_plot_Point_2D(z,'.k','-b',2,40);
sugr_plot_Point_2D(u,'.k','-b',2,40);
sugr_plot_Point_2D(v,'.k','-b',2,40);
sugr_plot_Line_2D(l,'-k','-r',2,40);
sugr_plot_Line_2D(m,'-k','-r',2,40);
axis equal
title('Intersection')

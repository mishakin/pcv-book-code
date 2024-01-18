%% DEMO 2D projective geometry
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 07/17
% wenzel@igg.uni-bonn.de

clc
close all
clearvars

addpath(genpath('../General-Functions'))

fprintf('\n------ DEMO 2D projective geometry ------\n')

fprintf('\nIn order to shorten the output, we display transposed vectors\n')
fprintf('\nGiven: 2D points\n')
x1 = [1;1;1]; disp(['x_1 = [', num2str(x1'),']'])
x2 = [1;2;1]; disp(['x_2 = [', num2str(x2'),']'])
x3 = [6;6;2]; disp(['x_3 = [', num2str(x3'),']'])

fprintf('\nlines l1 = x2 ^ x3 etc.\n')
l1 = calc_S(x2)*x3; disp(['l_1 = [', num2str(l1'),']'])
l2 = calc_S(x3)*x1; disp(['l_2 = [', num2str(l2'),']'])
l3 = calc_S(x1)*x2; disp(['l_3 = [', num2str(l3'),']'])

fprintf('\ny-axis\n')
ly = [1;0;0]; disp(['l_y = [', num2str(ly'),']'])

fprintf('\nintersection x4 = l3 ^ ly\n')
x4 = calc_S(ly)*l3;
disp(['x_4 = [', num2str(ly'),']  --> point at infinity in the direction of x-axis'])

plot_init;
figure('Color','w')
set(gca,'xlim',[0.5 3.5],'ylim',[0.5 3.5]);
hold on;

h = plot_2D_line(l1); set(h,'LineWidth',2);
h = plot_2D_line(l2); set(h,'LineWidth',2);
h = plot_2D_line(l3); set(h,'LineWidth',2);

x3 = x3/norm(x3);
scatter([x1(1),x2(1),x3(1)],[x1(2),x2(2),x3(2)],'ko','filled')


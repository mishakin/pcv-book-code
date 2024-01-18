%% eq_15_120_theoretical accuracy
%
% determines precision as a function of grid spacing
% plots grid and precision (10 levels)

% clear all
clearvars
close all
% type_data: type of generated data
% 15     4x4 Grid
%
% wf 7/2014

addpath(genpath('../../General-Functions'))
addpath('Functions')
addpath('Data')
plot_init;

disp('--- analyse theoretical precision ---')

init_rand = 6;
type_robust = 0
%             0 not robust
%             1 only dem
%             2 only points
%             3 both

out_C = 1;
out_print = 0;
Tol =0.15;

print_type = 0;
plot_type  = 0;

init_rand_seed(init_rand);

%% generate dem point cloud
sigma_z=zeros(10,1);
n=4
d_max=10;
for d=1:d_max
    [points,BB,delta_x,sigma_k,sigma_s,out_in,dem]=simulate_points_dem_15_flat(d,n);

    Np = size(points,1);


    % smooth non-robust

    type_robust=0;

    starttime = cputime
    [ds,S,Sigma,Np,Nr,Mc,ver,A,w,w_f,W] =...
        smooth_dem_robust_bilinear...
        (points,BB,delta_x,sigma_k,out_C,type_robust,out_in,...
        print_type,plot_type);
    complete_time_for_solution=cputime-starttime

    figure
    hold on
    mesh(S);
    title(strcat('Height = \sigma_z, d=',num2str(2*d)))
    view([-33,63])
    sigma_centre=sqrt(S(2+(n-1)*d,2+(n-1)*d));
    sigma_z(d)=sigma_centre;
end

std_z=sigma_z'

figure
hold on
d_range=1:d_max;
plot(2*(1:d_max),sigma_z(1:d_max),'ob','MarkerSize',12)
B=regress(sigma_z(d_range).^2,[ones(length(d_range),1),(2*d_range)'.^2])
plot(2*(1:d_max),sqrt(B(1)+B(2)*(2*(1:d_max))'.^2),'-r','LineWidth',2)
title(strcat('\sigma_z =',num2str(B(1)),' + ',num2str(B(2)),'d, (d > 1)'))
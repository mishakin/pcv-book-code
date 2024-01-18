%% Fig. 16.21 theoretical quality  of surface reconstruction
%
% Precision as a function of position
% Sensitivity wrt outliers in given heights
%
% wf 6/2015
%
% type_data: type of generated data
% 0     a set of given points


close all
% clear all
clearvars

addpath(genpath('../../General-Functions/'));
addpath('Functions')
addpath('Data')

ss = plot_init;



disp('---------------------------------------------------------')
disp('----- theoretical quality of surface reconstruction -----')
disp('---------------------------------------------------------')

%% set parameters
init_rand    =   2     % Example: with intit_rand = 2
N            =  25     % Example: with N = 25
sigma_h      =   1     % Example: with 1
sigma_k0     = 250     % Example: with 250
Resolution   =  50     % Example: with 50
sigma_k      = sigma_k0/Resolution   % adapted for resolution

out_sensitivity = 1  % = 0 if number of observations too large

type_robust = 0;
%             0 not robust
%             1 only dem
%             2 only points
%             3 both

out_C = 1;
out_print = 0;
Tol =0.15

print_type = 1;
plot_type  = 0;

init_rand_seed(init_rand);


%% generate dem point cloud

[points,BB,delta_x,dem,bi,dt]=simulate_points_dem_16_precision(N,sigma_h,Resolution);
out_in=ones(size(points,1),1);
Np = size(points,1);


%% interpolate

starttime = cputime;

[ds,S,Sigma,Np,Nr,Mc,ver,A,w,w_f,W] =...
    smooth_dem_robust_bilinear...
    (points,BB,delta_x,sigma_k,out_C,type_robust,out_in,...
    print_type,plot_type);

disp('---------------------------------------------------------')
complete_time_for_solution=cputime-starttime

%% sensitivity factors
rk=zeros(Np,1);
mu=zeros(Np,1);
for k=1:Np
    ak = A(k,:)';                          % k-th row of Jacobian
    sigma_k = ak'*Sigma*ak;                % standard deviation of k-th fitted height
    rk(k) = (1-sigma_k)/sigma_h^2;         % redundancy number
    mu(k) = sqrt(sigma_k/(1-sigma_k));     % sensitivity factor
end


%%
X=([0:Nr-1]'*ones(1,Mc))*delta_x+BB(1);
Y=(ones(Nr,1)*[0:Mc-1])*delta_x+BB(2);
%% plot fitted dem

figure('Color','w','Position',[20,0.2*ss(2),0.6*ss(1),0.6*ss(2)])
hold on
%subplot(2,2,1)

hold on
surf(X,Y,ds); %,'EdgeColor',[0,0,0])% ,129*ones(size(S))
colormap(gray)

title('fitted dem: z')
alpha(0.3)
col=zeros(Np,2);
for n=1:Np
    plot3(points(n,1),points(n,2),points(n,3),'.k','MarkerSize',15)
end
view([-29,65])
xlim([BB(1),BB(3)])
ylim([BB(2),BB(4)])
zlim([0,10])
axis equal

%% plot theoretical quality
figure('Color','w','Position',[30,0.15*ss(2),0.6*ss(1),0.6*ss(2)])
hold on
Smax=max(sqrt(S(:)));
surf(X,Y,sqrt(S)/Smax*10,'EdgeColor',[0,0,0])
colormap(cool);
title('standard deviations $\sigma_z$ - sensitifyty factors $\mu$')
col=zeros(Np,2);
alpha(0.3)
z=zeros(Np,1);
for n=1:Np
    z(n)=interpolate_bilinear(S,points(n,1),points(n,2),delta_x,BB,sigma_h);
    plot3(points(n,1),points(n,2),sqrt(z(n)),'.k','MarkerSize',15)
end
if out_sensitivity ~= 0
    for n=1:Np
        plot3([points(n,1),points(n,1)],[points(n,2),points(n,2)],[0,mu(n)],'-r','LineWidth',2)
    end
end
plot3(0,0,0,'.k')
view([-29,65])
xlim([BB(1),BB(3)])
ylim([BB(2),BB(4)])
zlim([-1,7])
axis equal

x_y_sigma_std_fitted_points = [points(:,[1,2,4]),sqrt(z)];

Std_z=sqrt(S);

%% plot sensitivity factors

if out_sensitivity ~= 0
    figure('Color','w','Position',[40,0.20*ss(2),0.5*ss(1),0.5*ss(2)])
    hold on
    R=Resolution;
    plot([0,R+5,R+5,0,0],[0,0,R+5,R+5,0],'-k')
    for k=1:Np
        i=points(k,1)+1;
        j=points(k,2)+1;
        plot(i,j,'.b','MarkerSize',10)
        text(i+0.5,j+0.5,num2str(round(mu(k)*10.001)/10,'%6.1f'))
    end
    disp('---------------------------------------------------------')
    title('sensitivity factors')
    axis equal
    axis off

    disp('             sensitivity                ')
    disp('         i         j        rk        mu')
    [points(:,1:2),rk,mu]
end
average_redundancy_observations=sum(rk)/Np


return

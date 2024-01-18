%% Perturb sampled 3D 2D points pairs
%
% [X,y] = sugr_perturb_3D_2D_point_pairs(Xt,yt,sigma_X,sigma_y,fp);
% 
% Xt,yt     = true point pairs, sugr-objects, 3D and 2D points, resp.
% sigma_x   = std for perturbation of X's: CovM = sigma_x^2 I_3
% sigma_y   = std for perturbation of y's: CovM = sigma_y^2 I_2
% fp        = magnification factor for plot
%
% [X,y]      perturbed point pairs, sugr-objects, 3D and 2D points, resp.
%
% Wolfgang Förstner 02/2013
% wfoerstn@uni-bonn.de
%
% See also sugr_perturb_2D_point_pairs, sugr_perturb_2D_point_pairs_spherical, 
% sugr_perturb_Lines_2D

function [X,y] = sugr_perturb_3D_2D_point_pairs(Xt,yt,sigma_X,sigma_y,fp)

global plot_option 

Xte   = Xt.e;
CXXe  = sigma_X^2*eye(3);
yte   = yt.e;
Cyye  = sigma_y^2*eye(2);

N = size(Xte,1);

X.e   = zeros(N,3);
X.Cee = zeros(N,3,3);
y.e   = zeros(N,2);
y.Cee = zeros(N,2,2);
for n=1:N
    % true Euclidean coordinates
    X_true_e = Xte(n,1:3)';
    y_true_e = yte(n,1:2)';
    % perturbed point pair
    X.e(n,:)       = sugr_rand_gauss(X_true_e,CXXe, 1)';
    y.e(n,:)       = sugr_rand_gauss(y_true_e,Cyye, 1)';
    X.Cee(n,:,:) = CXXe;
    y.Cee(n,:,:) = Cyye;
end
X.e;
y.e;

% plot
ss = plot_init;
if plot_option > 0
    figure('Color','w','Position',[ss(1)/4,ss(2)/3,ss(1)/2,ss(2)/2])
   
    subplot(1,2,1)
    hold on
    for n=1:N
        plot3(X.e(n,1),X.e(n,2),X.e(n,3),'.b');
    end
    title('3D points')
    axis equal
    
    subplot(1,2,2)
    hold on
    % 2D 
    for n=1:N
        ye = y.e(n,:)';
        Cyy = squeeze(y.Cee(n,:,:));
        yn = sugr_Point_2D(ye,Cyy);
        sugr_plot_Point_2D(yn,'.k','-b',2,fp);
    end
    title('image points')
    axis equal
end

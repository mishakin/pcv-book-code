%% Perturb sampled 2D points pairs sphericaly
%
% PP = sugr_perturb_2D_point_pairs_spherical(PP_true,sigma_x,sigma_y,rho,fp);
% 
% PP_true = struct, sugr point pairs 2D, true point pairs
%              .h(N,6)  homogeneous coordinates of N point pairs stored in each row [x1_n', x2_n']
%              .Crr(N,4,4) according reduced CovM
% sigma_x   = standard deviation for perturbation of x's: CovM = sigma_x^2 I_2
% sigma_y   = standard deviation for perturbation of y's: CovM = sigma_x^2 I_2
% rho       = correlation coefficient
% fp        = magnification factor for plot
%
% PP perturbed point pairs
%
% Wolfgang Förstner 10/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_perturb_2D_point_pairs, sugr_perturb_3D_2D_point_pairs,
% sugr_perturb_Lines_2D

function PP = sugr_perturb_2D_point_pairs_spherical(PP_true,sigma_x,sigma_y,rho,fp)

global plot_option 

Pth   = PP_true.h;
% PtCrr = PP_true.Crr;

N = size(Pth,1);

Cppr     = [sigma_x^2*eye(2) sigma_x*sigma_y*rho*eye(2);...
            sigma_x*sigma_y*rho*eye(2) sigma_y^2*eye(2)];
for n=1:N
    % perturbed point pair
    dp       = sugr_rand_gauss(zeros(4,1),Cppr, 1);
    ps(1:3)  = sugr_ghm_update_vector(Pth(n,1:3)',dp(1:2));
    ps(4:6)  = sugr_ghm_update_vector(Pth(n,4:6)',dp(3:4));
    PP.h(n,:)     = ps;
    PP.Crr(n,:,:) = Cppr;
    PP.type(n)    = 8;
end

if plot_option > 0
    figure
    for n=1:N
         xh = PP.h(n,1:3)';
         yh = PP.h(n,4:6)';
         Cxx = null(xh') * squeeze(PP.Crr(n,1:2,1:2)) * null(xh')';
         Cyy = null(yh') * squeeze(PP.Crr(n,3:4,3:4)) * null(yh')';
         x = sugr_Point_2D(xh,Cxx);
         y = sugr_Point_2D(yh,Cyy);
         sugr_plot_Point_2D(x,'.k','-r',2,fp);
         sugr_plot_Point_2D(y,'.k','-b',2,fp);
         l=sugr_construct_join_Line_2D(x,y);
         sugr_plot_Line_2D(l,'-k','-w',1,5,[0,0]);
    end
    axis equal
end


%% Perturb sampled 2D points pairs
%
% PP = sugr_perturb_2D_point_pairs(PP_true,sigma_x,sigma_y,rho,fp);
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
% Wolfgang Förstner 2/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_perturb_2D_point_pairs_spherical, sugr_perturb_3D_2D_point_pairs,
% sugr_perturb_Lines_2D

function PP = sugr_perturb_2D_point_pairs(PP_true,sigma_x,sigma_y,rho,fp)

global plot_option

Pth   = PP_true.h;

N = size(Pth,1);

Cppe     = [sigma_x^2*eye(2) sigma_x*sigma_y*rho*eye(2);...
    sigma_x*sigma_y*rho*eye(2) sigma_y^2*eye(2)];
for n=1:N
    % true Euclidean coordinates
    x_true_e = Pth(n,1:2)'/Pth(n,3);
    y_true_e = Pth(n,4:5)'/Pth(n,6);
    p_true_e = [x_true_e;y_true_e];
    % perturbed point pair
    pe       = sugr_rand_gauss(p_true_e,Cppe, 1);
    %pp       = sugr_Point_Pair_2D(pe,Cppe);
    %     PP.h(n,:)     = pp.h';
    %     PP.Crr(n,:,:) = pp.Crr;
    %     PP.type(n)    = 8;
    %
    % homogeneous coordinates
    x1h = [pe(1:2);1];
    n1  = norm(x1h);
    x2h = [pe(3:4);1];
    n2  = norm(x2h);
    ph   = [x1h/n1;x2h/n2];
    % covariance matrix
    Chh  = [Cppe(1:2,1:2) zeros(2,1) Cppe(1:2,3:4) zeros(2,1);...
        zeros(1,6);...
        Cppe(3:4,1:2) zeros(2,1) Cppe(3:4,3:4) zeros(2,1);...
        zeros(1,6)];
    J    = [null(x1h')/n1, zeros(3,2); ...
        zeros(3,2), null(x2h')/n2];
    Crr    = J' * Chh * J;               % CovM of reduced vector
    %
    PP.h(n,:)     = ph';
    PP.Crr(n,:,:) = Crr;
    
end


for n=1:N
    xh = PP.h(n,1:3)';
    yh = PP.h(n,4:6)';
    Cxx = null(xh') * squeeze(PP.Crr(n,1:2,1:2)) * null(xh')';
    Cyy = null(yh') * squeeze(PP.Crr(n,3:4,3:4)) * null(yh')';
    x = sugr_Point_2D(xh,Cxx);
    y = sugr_Point_2D(yh,Cyy);
    if plot_option > 0
        sugr_plot_Point_2D(x,'.k','-r',2,fp);
        sugr_plot_Point_2D(y,'.k','-b',2,fp);
        l = sugr_construct_join_Line_2D(x,y);
        sugr_plot_Line_2D(l,'-k','-w',1,5,[0,0]);
    end
    
end



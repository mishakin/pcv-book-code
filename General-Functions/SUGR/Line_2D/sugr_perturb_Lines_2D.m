%% Perturbes given image lines through one point
%
% lstruct = sugr_perturb_Lines_2D(x0,d,lines,resolution,sigma, rigorous)
%
% * xo = [x,y,phi] true line centroid, true direction of line
% * d   = length of line
% * lines = Nx4 [u0,v0,phi_0,lengthl] of N true lines
% * resolution = sampling distance for ploting
% * sigma [pel] = standard deviation of single pixels
%      sigma_q = sigma_0/length   sigma of cross-position deviation
%      sigma_phi^2 = sigma_0^2*length^3/12
% * rigorous = 0 approximate covariance matrix, = 1 rogorous covariance matrix
%
% lstruct = lstruct.h    =  N x 3 matrix of homogeneous coordinates
%           lstruct.Crr  =  N x 2 x 2 matrix of reduced covariance matrices
%           lstruct.type =  N x 1 vector of type = 2        
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
% 
% See also sugr_Line_2D, sugr_generate_true_Lines_2D_one_Point

function lstruct = sugr_perturb_Lines_2D(x0,d,lines,resolution,sigma, rigorous) %#ok<INUSD>


global plot_option
global print_option

%% Initiate

[N,~] = size(lines);
c = 1;

if plot_option > 0
    van_1 = x0/norm(x0);
    van_1(3)=van_1(3)*c;
    %plot(2*d*[1,-1,-1,1,1]',2*d*[1,1,-1,-1,1]','-b');
    xlim([-2.1,2.1]);
    ylim([-2,2.7])
    if abs(van_1(3)) > sqrt(1/2)*d* norm(van_1(1:2))
        van_1_e=van_1(1:2)/van_1(3);
        plot(van_1_e(1),van_1_e(2),'or')
    end
end

%% prepare generation of lines

lstruct.h    = zeros(N,3);
lstruct.Crr  = zeros(N,2,2);
lstruct.type = 2 * ones(N,1);

ii = 0;

Delta_x=d/resolution;
sigma=sigma*Delta_x;

% colors
% col(0+1)='g';
% col(1+1)='r';

% maxl=0;
% minl=10000;

%% generate N lines
for n = 1:N
    % generate line to point van_1

    % Generate noisy data
    % centre
    z0     = lines(n,1)+1i*lines(n,2);
    
    % true direction
    phi_0  = lines(n,3);

    % length
    lengthl= lines(n,4);
    
    N_points = round(lengthl);

    % std. deviation across line
    sigma_k = sigma/sqrt(N_points);
    
    % noisy centre
    z0_n = z0 + (randn(1)+1i*randn(1)) * sigma_k;
    
    % std. dev. direction
    sigma_phi= sigma/sqrt(N_points^3/12)/Delta_x;
    
    % noisy direction 
    phi_n = phi_0 + randn(1)*sigma_phi;

    % start and end point
    dz_n = lengthl/2*exp(1i*phi_n)/resolution*d;

    zs_0 = z0_n - dz_n;
    ze_0 = z0_n + dz_n;

    % convert to 2D
    xs = real(zs_0) ;
    ys = imag(zs_0) ;
    xe = real(ze_0) ;
    ye = imag(ze_0) ;
%     [xs,ys,xe,ye,phi_n];
    ii=ii+1;
    Line                = sugr_Line_2D((xs+xe)/2,(ys+ye)/2,phi_n-pi/2,sigma_phi,sigma_k);
    lstruct.h(ii,:)     = Line.h';
    lstruct.Crr(ii,:,:) = Line.Crr;
end;

% [maxl,minl];
nii=ii;

%% Check and Plot if plot_option fulfilled

if print_option > 0
    Omega_d = 0;
    Omega_p = 0;
    for j=1:nii
        l       = sugr_select_Line_2D(lstruct,j);
        Clhh    = sugr_get_CovM_homogeneous_Vector(l);
        var_c   = van_1' * Clhh * van_1;
        c       = l.h' * van_1;
        Omega_d = Omega_d + c^2;
        Omega_p = Omega_p + c^2 / var_c;
    end;
    % check
%     RMS_d=sqrt(Omega_d/nii)*resolution;
%     sigma_0_est = sqrt(Omega_p/(nii-2));
end

if plot_option > 0
    for j=1:nii        
        l       = sugr_select_Line_2D(lstruct,j);
        sugr_plot_Line_2D(l,'-w','-k',2,100,[0,0]);
    end
    axis equal
end



            
  





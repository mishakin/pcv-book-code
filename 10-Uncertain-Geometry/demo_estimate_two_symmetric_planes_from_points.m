%% test_two_symmetric_planes_from_points
%
% example PCV Sect. 10.6.4, Fig. 10.29
%
% The parameter model_II is used to choose the estimation model
%
% % Wolfgang Förstner 7/2012
% wfoerstn@uni-bonn.de

%% 
close all
clearvars
clc

addpath(genpath('../General-Functions'))
addpath(genpath('Functions'))

global min_redundancy
global plot_option
global print_option_estimation
ss = plot_init;

disp('---------------------------------------')
disp(' Estimate parameters of symmetric roof ')
disp(' example PCV Sect. 10.6.4, Fig. 10.29  ')
disp('---------------------------------------')

%% Choose estimation model

% model_II = 0; % One step estimation GHM with 
model_II = 1; % Estimate individually, Model D + C

switch model_II
    case 0
        disp('One step estimation GHM with constraints')
    case 1
        disp('Estimate individually, Model D + C')
end

% read coordinates according to PCV?
read_PCV_coordinates = 0; % coordinates from book, only 4-digits
% read_PCV_coordinates = 0; % generate random points

%% initialize 
factor = 2;                  % for scaling the complete data set
slope = -0.6;                % slope of roofs
NL    = 5;                   % number of points on left roof
NR    = 6;                   % number of points on left roof
sigma =  factor*0.05;        % std. dev. of points
maxiter = 5;                 % iteration for constraints
Tol     = 0.0001;
disp(strcat('Tolerance for convergence :', num2str(Tol)))

rangeL = factor*[-1,0,0,2]; % x-range, y-range of roofs
rangeR = factor*[ 0,1,0,2];

plot_option = 1;
print_option_estimation = 2;

% initialize SUGR
sugr_INIT
min_redundancy = 100000;

init_rand      = 3;          % random seed (default = 3)
init_rand = init_rand_seed(init_rand);
%% Generate planes
disp('Generate planes...');
% common horizontal line: height 1.5, in Y-direction
X = [0,0,factor*1.5,1]';             % start points
Y = [0,factor*1,factor*1.5,1]';      % end point
L = calc_Pi(X)*Y;                    % roof line
L_true = L;
disp(['   Roof edge L = [',num2str(L'),']'])

% slope, tan(alpha) 
X_L = X + factor*[1,0,-slope,0]';
X_R = X + factor*[-1,0,-slope,0]';
AL  = +calc_Pidual(X_L)'*L;
AR  = -calc_Pidual(X_R)'*L;
AL_true = AL;
AR_true = AR;
disp(['   Left roof plane AL  : [',num2str(AL'),']'])
disp(['   Right roof plane AR : [',num2str(AR'),']'])

%% generate points
disp(' ')
disp('Generate points...');
if plot_option > 0    
    figure('Color','w','Position',[0.2*ss(1), 0.35*ss(2),ss(1)/2,ss(2)/2])
    hold on
    axis equal
end

disp('---Generate points on left plane')
XL = generate_points_on_plane(NL,AL,sigma,rangeL,0,read_PCV_coordinates);
disp('---Generate points on right plane')
XR = generate_points_on_plane(NR,AR,sigma,rangeR,1,read_PCV_coordinates);

xlim(factor*[-1.5,1.5]);
   ylim(factor*[-0.5,2.5]);
   zlim(factor*[0,2.5]);
   axis equal
set(gca,'CameraPosition',[-14,31,27])

%%  approximate values of the two planes: algebraic solution X_i' * A = 0
[~,~,V] = svd(XL.h);  % 
ALa  =V(:,4);
[U,D,V] = svd(XR.h);
ARa = V(:,4);

planes_appr = [ALa'/ALa(3);ARa'/ARa(3)];
disp(['Approximate plane parameters left  : [',num2str(planes_appr(1,:)),']'])
disp(['Approximate plane parameters right : [',num2str(planes_appr(2,:)),']'])

% Possibly show coordinates
% XL.h(:,1:3)./(XL.h(:,4)*[1,1,1]);
% XR.h(:,1:3)./(XR.h(:,4)*[1,1,1]);

%% Estimate complete roof
disp(' ')
disp('Perform estimation...');
if model_II==1
    
    disp('-------- Estimate individually, Model D + C');
    
    %% First step: estimate individually using model D
    disp(' ')
    disp('ML-estimation left plane ---------------------------------------')
    [AL_est,sigma_0L,RL] = sugr_estimation_ml_Plane_from_points(XL,ALa,0.01,5);
    disp(' ')
    disp('ML-estimation right plane --------------------------------------')
    [AR_est,sigma_0R,RR] = sugr_estimation_ml_Plane_from_points(XR,ARa,0.01,5);
    disp(['   left plane: R = ',num2str(RL), ', sigma_0 = ',num2str(sigma_0L)]);
    disp(['   right plane: R = ',num2str(RR), ', sigma_0 = ',num2str(sigma_0R)]);
    
%     AL_est.h/AL_est.h(3)
%     AR_est.h/AR_est.h(3)
    
    disp(['   Estimated left plane  (1): ',num2str(AL_est.h'/AL_est.h(3))])
    disp(['   Estimated right plane (1): ',num2str(AR_est.h'/AR_est.h(3))])
    
    Crr  = [AL_est.Crr zeros(3)   ;...
            zeros(3)   AR_est.Crr];
%     Crr_x_10000 = Crr*10000 
    
    %% Second step: Impose constraints on plane parameters using model C
    % see PCV Sect. 4.8.2.6 
    disp('-------- -------- Impose constraints: Model C');
    Nm  = zeros(2);
    hv  = zeros(2,1);
    l  = [AL_est.h;AR_est.h];   % observations == plane parameters
    l0 = l;
    AL = l(1:4);
    AR  =l(5:8); 
    for iter = 1:maxiter

        disp(['-------- -------- -------- C Iteration ',num2str(iter),'----------------'])

        AL0 = l0(1:4);                    % left plane A
        AR0 = l0(5:8);                    % right plane B
        est_L0 = calc_Pidual(AL0)*AR0;  % roof edge L
        eL0C = est_L0'/norm(est_L0(1:3));

        % Jacobian of constraints, see PCV (10.318)
        Bt=[2*[-AL0(1)*AR0(3)^2,...
            -AL0(2)*AR0(3)^2,...
            +AL0(3)*(AR0(1)^2+AR0(2)^2),...
            0]*null(AL0'),...
            2*[+AR0(1)*AL0(3)^2,...
            +AR0(2)*AL0(3)^2,...
            -AR0(3)*(AL0(1)^2+AL0(2)^2),...
            0]*null(AR0'),
            [[+AR0(2),-AR0(1),0,0]*null(AL0'),...
            [-AL0(2),+AL0(1),0,0]*null(AR0')]];                            %#ok<COMNL>
        % Jacobian for reducing the plane coordinates
        Jrh = [null(l0(1:4)')' zeros(3,4); zeros(3,4) null(l0(5:8)')'];
        % constraints,  see PCV (10.318) 
        cg = [-(AL0(3)^2*(AR0(1)^2+AR0(2)^2)-AR0(3)^2*(AL0(1)^2+AL0(2)^2));...
            -(AL0(1)*AR0(2)-AR0(1)*AL0(2))] + Bt*Jrh*(l0-l);
        % normal equation matrix (4.463)
        Nm = Bt*Crr*Bt';
        % correction to fitted observations (4.465)
        Dlr = Crr * Bt' * inv(Nm) * cg + ...
            Jrh * l;                                                       %#ok<MINV>
        disp(['Correction of observations (reduced plane parameters):', num2str(Dlr')])
        
        % corrected fitted observations,  see PCV (10.255)
        l0(1:4)=sugr_ghm_update_vector(l0(1:4),Dlr(1:3)); 
        l0(5:8)=sugr_ghm_update_vector(l0(5:8),Dlr(4:6));
        
        if (abs(Dlr(:))./sqrt(diag(Crr)) < Tol) | (iter==maxiter)         %#ok<OR2,BDSCA>
            disp('-------- -------- -------- after last iteration')
            % covariance matrix of fitted observations
            % CovM(est_l) = CovM(l) - CovM(v), see (4.68) with (4.467)
            estCrr = Crr - Crr * Bt' * inv(Nm) * Bt * Crr;                 %#ok<MINV>
            ver = Jrh * (l0-l);
            Omega    = ver' * inv(Crr) * ver;                              %#ok<MINV>
            sigma_0 = sqrt(Omega/2);
            % redundancy is R=2 (2 constraints)
            disp(['    R = 2 sigma_0 = ',num2str(sigma_0)]);
            break
        end
    end
    est_AL = l0(1:4);
    est_AR = l0(5:8);

    disp(' ')
    disp(['Estimated left plane  : [',num2str(est_AL'/est_AL(3)),']'])
    disp(['Estimated lright plane: [',num2str(est_AR'/est_AR(3)),']'])

    est_L = calc_Pidual(est_AL)*est_AR;
    
    disp(['Estimated roof line   : [',num2str(est_L'/norm(est_L(1:3))),']'])

    % covariance matrix of planes (from reduced to homogeneous coord.)
    CestAestA = Jrh' * estCrr * Jrh;

    % covariance matrix of 3D line
    JLA = [-calc_Pidual(est_AR),calc_Pidual(est_AL)];
    CLL = JLA * CestAestA * JLA';
        
    % take vector L(1:3) as point in P^2, spherically normalized
    dir = sugr_Point_2D(est_L(1:3),CLL(1:3,1:3));
    % eigenvalues yield uncertainty of direction vector, 
    % one ev is 0, the slope of the fitter roof edge
    % the other is the variance of the azimuth
    ev = eigs(dir.Crr)*180/pi;
    disp(['Std.dev. of azimuth (degrees): ', num2str(sqrt(max(ev(:))))])
    
else
    % one step estimation GHM with constraints
    disp('-------- One-step estimation Model E')
    
    [estA,estX,sigma_0,R] = ...
        sugr_estimation_ml_symmetric_Planes_from_points...
        (XL,XR,[ALa;ARa],Tol,maxiter);

    disp('-------- -------- -------- after last iteration')
    
    disp(['Estimated sigma_0^2  :',num2str(sigma_0),'^2'])
    % fitted planes and points
    AL0 = estA.h(1:4);
    AR0 = estA.h(5:8);
    XL0 = estX(1:NL,:);
    XR0 = estX(NL+1:NL+NR,:);
    est_AL = AL0;
    est_AR = AR0;
    
    XLe = XL0(:,1:4)./(XL0(:,4)*[1,1,1,1]);
    XRe = XR0(:,1:4)./(XR0(:,4)*[1,1,1,1]);
%     XLe*AL0/norm(AL0(1:3));
%     XRe*AR0/norm(AR0(1:3));

    eAL = est_AL'/est_AL(3);
    eAR = est_AR'/est_AR(3);
    disp(['Estimated left plane  :',num2str(eAL)])
    disp(['Estimated right plane :',num2str(eAR)])

    est_L = calc_Pidual(est_AL)*est_AR;
    estL = est_L'/norm(est_L(1:3));
    disp(['Estimated roof line   :',num2str(estL)])

    % covariance matrix of 3D line
    JLA = [-calc_Pidual(est_AR),calc_Pidual(est_AL)];
    CLL = JLA * estA.Cxx * JLA';
    
    % take elements L(1:3) as direction on S^2 (see above)
    dir = sugr_Point_2D(est_L(1:3),CLL(1:3,1:3));
    ev = eigs(dir.Crr)*180/pi;
    disp(['Std.dev. of azimuth (degrees): ', num2str(sqrt(max(ev(:))))])
    
%%    
if plot_option > 0
    figure('Color','w','Position',[0.2*ss(1), 0.1*ss(2),ss(1)/2,ss(2)/2])
    hold on
    r = rangeL;
    A = AL_true;
    z = -(A(1)*r(1:2) + A(2)*r(3:4) + A(4)) / A(3);
   
   plot3([r(1),r(2),r(2),r(1),r(1)],[r(3),r(3),r(4),r(4),r(3)],[z(1),z(2),z(2),z(1),z(1)],'--k');
   for n=1:NL
       plot3(XLe(n,1),XLe(n,2),XLe(n,3),'.k','MarkerSize',15);
   end
   
   r=rangeR;
   A=AR_true;
   z=-(A(1)*r(1:2) + A(2)*r(3:4) + A(4)) / A(3);
   plot3([r(1),r(2),r(2),r(1),r(1)],[r(3),r(3),r(4),r(4),r(3)],[z(1),z(2),z(2),z(1),z(1)],'--k');
   
   for n=1:NR
       plot3(XRe(n,1),XRe(n,2),XRe(n,3),'.k','MarkerSize',15);
   end
  
   xlim(factor*[-1.3,1.3]);
   ylim(factor*[-0.3,2.3]);
   zlim(factor*[0,2.5]);
   xlabel('X')
   ylabel('Y')
   zlabel('Z')
   axis equal
   
   set(gca,'CameraPosition',[-14,31,27])
end
   
end

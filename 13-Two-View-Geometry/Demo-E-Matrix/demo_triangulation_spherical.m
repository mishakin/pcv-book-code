%% demo triangulation
%
% simulate:
% X        = 4x1 vector
% [b,R]    = relative orientation
% sigma    = radial error odf direction
% max_iter = maximum iterations
% Tol      = tolerance for stopping iterations
% k        = critical value
%
% derive 
% u, v     = 3x1 vectors, spherically normalized
%
% evaluate estimated 3D point with its covaraince matrix
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de


% clear all
clearvars
close all

addpath(genpath('../../General-Functions'))
addpath('../Functions')

sugr_INIT

ss = plot_init();

disp(' --------------------------------------');
disp(' ------ Triangulation (Alg. 21) -------');
disp(' --------------------------------------');

N_samples = 1000;

sigma = 0.0001;

% random seed
% init_rand      = 1; 
init_rand      = 215741;

disp(['Number of samples                      : ', num2str(N_samples)])
disp(['Standard deviation of directions [rad] : ', num2str(sigma)])

%% random number intialization
init_randn = init_rand_seed(init_rand);
disp(['Random seed                            : ' , num2str(init_randn)])

%%   
fcond = 1;
X_true = [[0.5,1.0,1.5]*fcond,1]';
Jr_true = null(X_true');
% b = [1,0,0]'; 
% b = [0,0,1]'; 
b = randn(3,1);
b = b/norm(b);

xi = X_true(4)/norm(X_true(1:3));
infinite = xi < 0.0001; % only for plotting

Rot = calc_Rot_r(randn(3,1));
sigma = 0.0001;
max_iter = 10;
Tol = 0.0001;
k = 10;
print_option = N_samples == 1; 

X_est_all = zeros(N_samples,3);
C_est_all = zeros(N_samples,3,3);
s_all     = zeros(N_samples,1);

%% for all samples
t = cputime;
for n = 1:N_samples
    
    % derive u, v
    [u,v] = generate_point_pair(X_true,b,Rot,sigma);
    
    % E        = sugr_SkewMatrix(b)*Rot';
    % sigma_c  = sigma*sqrt(u'*E*E'*u+v'*E'*E*v);
    % cgl      = u'*E*v;  
    % s_all(n) = abs(cgl)/sigma_c;
    % estimate X
    
    [X_est, est_sigma_0, est_u, est_v, f] = ...
        sugr_triangulate_two_spherical_cameras(...
        b, Rot, u, v, sigma, Tol, max_iter, k, print_option);
    
    X_est_all(n,:) = (Jr_true'*X_est.h)';
    
    s_all(n)       = est_sigma_0; 
end

%%
disp(['CPU time for ', num2str(N_samples),' samples              : ',num2str(cputime-t)]);
m_X      = mean(X_est_all);    
cov_X    = cov(X_est_all);
disp(['mean reduced coordinates :', num2str(m_X)]);

%% Plot
% plotting of not too far
if ~infinite
    figure('name','est.cov','Color','w', 'Position',[100 0.52*ss(2) ss(1)/3 0.4*ss(2)])
    hold on
    plot3(X_true(1),X_true(2),X_true(3),'or');
    plot3([0,b(1)],[0,b(2)],[0,b(3)],'-b','LineWidth',3);
    plot3([0,X_true(1)],[0,X_true(2)],[0,X_true(3)],'--r','LineWidth',2);
    plot3([b(1),X_true(1)],[b(2),X_true(2)],[b(3),X_true(3)],'--r','LineWidth',2);
    box on
    title('estimated covariance matrix')
    axis equal

    Xemp.h = X_est.h;
    Xemp.Crr = cov_X;
    [~,Cee] = sugr_get_Euclidean_Point_3D(Xemp);
    [R,D] = eig(Cee);

    f = 0.4/sqrt(D(1,1));
    [X,Y,Z] = ellipsoid(0,0,0,sqrt(D(1,1))*f,sqrt(D(2,2))'*f,sqrt(D(3,3))*f);
    K = [X(:),Y(:),Z(:)]*R'+ones(441,1)*X_true(1:3)';

    surfl(reshape(K(:,1),21,21), reshape(K(:,2),21,21), reshape(K(:,3),21,21))
    colormap cool %copper
    axis equal

    figure('name','true.cov','Color','w', 'Position',[300+ss(1)/3 0.52*ss(2) ss(1)/3 0.4*ss(2)])
    hold on
    plot3(X_true(1),X_true(2),X_true(3),'or');
    plot3([0,b(1)],[0,b(2)],[0,b(3)],'-b','LineWidth',3);
    plot3([0,X_true(1)],[0,X_true(2)],[0,X_true(3)],'--r','LineWidth',2);
    plot3([b(1),X_true(1)],[b(2),X_true(2)],[b(3),X_true(3)],'--r','LineWidth',2);
    box on
    title('true covariance matrix')
    axis equal
        
    [~,Cee] = sugr_get_Euclidean_Point_3D(X_est);
    [R,D] = eig(Cee);

    f = 0.4/sqrt(D(1,1));

    [X,Y,Z] = ellipsoid(0,0,0,sqrt(D(1,1))*f,sqrt(D(2,2))'*f,sqrt(D(3,3))*f);
    K = [X(:),Y(:),Z(:)]*R'+ones(441,1)*X_true(1:3)';

    surfl(reshape(K(:,1),21,21), reshape(K(:,2),21,21), reshape(K(:,3),21,21))
    colormap copper
    axis equal

end

%% Analyse
check_estimation_result(1,zeros(3,1),X_est.Crr,s_all.^2,X_est_all,0.999,'Triangulation');
set(gcf,'Position',[100 20 ss(1)/2 0.4*ss(2)])


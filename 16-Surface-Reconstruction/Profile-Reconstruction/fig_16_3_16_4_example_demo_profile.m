%% Fig. 16.3 page 731 and 16.4 page 732
% demo profile reconstruction with differente priors.
%
% Wolfgang Förstner 2014-10-05
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

addpath(genpath('../../General-Functions/'));
addpath('Functions');

close all

%% plot settings
ss = plot_init;

%% initialize random number generation by fixed seed
init_rand_seed(4);

disp(' --------------------------------------------------')
disp(' ------ Fig. 16.3 demo profile reconstruction -----')
disp(' --------------------------------------------------')


%% generate true profile

disp(' ------ generate true profile -----')
N = 200;   % number of points total
N1 = 70;   % number of points in first (flat) part of the profil
N2 = N-N1; % remaing points

% first AR(1) - flat part of the profil
sigma_e = 1.0;  % process noise
sigma_n = 0.5;  % observation noise
a = 0.98;       % decay

x1 = zeros(N1,1);   % true profile
y1 = zeros(N1,1);   % observed sample points
for n=1+1:N1
    x1(n) = sum(a.*x1(n-1:-1:n-1))+randn(1)*sigma_e;
    % discontinuity at position 50
    if n==50
        x1(n)=x1(n)+8;
    end
    y1(n) = x1(n)+randn(1)*sigma_n;
end
% enforce slope 0
for i=1:N1
    y1(i)=y1(i) -  (i-1)*(x1(N1)-x1(1))/(N1-1);
    x1(i)=x1(i) -  (i-1)*(x1(N1)-x1(1))/(N1-1);
end

% second AR(2) - valley
sigma_e = 0.5;  % process noise
sigma_n = 0.5;  % observation noise
a = [2,-1]';

x2=zeros(N2,1); % true profile
y2=zeros(N2,1); % observed sample points
for n=2+1:N2
    x2(n) = sum(a.*x2(n-1:-1:n-2))+randn(1)*sigma_e;
    % discontinuity at position 50
    if n==50
        x2(n)=x2(n)+5.0;
    end
    y2(n) = x2(n)+randn(1)*sigma_n;
end
% enforce slope 0
for i=1:N2
    y2(i)=y2(i) -  (i-1)*(x2(N2)-x2(1))/(N2-1);
    x2(i)=x2(i) -  (i-1)*(x2(N2)-x2(1))/(N2-1);
end
% complete true profile
x = [x1;x2+x1(N1)];
% complete set of observed sample points
y = [y1;y2+y1(N1)];

% plot the profile
figure('name','Fig. 16.3: Profile Reconstruction','color','w',...
    'Position',[0.1*ss(1),0.1*ss(2),0.8*ss(1),0.8*ss(2)]);

ymin = min(y)-20; ymax = max(y)+20;
subplot(4,3,1); hold on;
plot(1:N,x,'-k','LineWidth',2)
ylim([ymin,ymax])
title('true profile')

%% sample
disp(' ------ sample profile -----')
Mset = [6,20,48,60,100,144,164,182]';
M = length(Mset);
select = zeros(8,1);
xs = zeros(8,1); ys = zeros(8,1);
for m = 1:M
    select(m) = Mset(m);
    xs(m) = x(Mset(m));
    ys(m) = y(Mset(m));
end

subplot(4,3,2);hold on
plot(1:N,x,'-k','LineWidth',2)
plot(Mset,ys,'.b','MarkerSize',15)
ylim([ymin,ymax])
title('true profile with sampled points')

subplot(4,3,3); hold on
plot(Mset,ys,'.b','MarkerSize',15)
ylim([ymin,ymax])
title('sampled points')

%% interpolation polynomial 0

disp(' ------ interpolation polynomial 0 -----')
yest = mean(ys);

subplot(4,3,4); hold on
plot(1:N,yest*ones(N,1),'-r','LineWidth',2)
plot(Mset,ys,'.b','MarkerSize',15)
plot(1:N,x,'--k','LineWidth',1)
ylim([ymin,ymax])
title('best horizontal line')

%% interpolation polynomial of order 3

disp(' ------ interpolation order-3-polynomial -----')
order = 3;
% coefficients of sampled coordinates
Xmatrix = [ones(M,1),Mset,Mset.^2,Mset.^3];
% estimate params by regression
estpar  = regress(ys,Xmatrix(:,1:order+1));

% reconstruction
xc = (1:N)';
Xpred=[ones(N,1),xc,xc.^2,xc.^3];
xpred = Xpred(:,1:order+1)*estpar;

subplot(4,3,5); hold on
plot(1:N,xpred,'-r','LineWidth',2)   % reconstructed polynomial
plot(Mset,ys,'.b','MarkerSize',15)   % sample points
plot(1:N,x,'--k','LineWidth',1)      % true profile
ylim([ymin,ymax])
title(strcat('best order-',num2str(order),'-polynomial'))

%% interpolation polynomial of order 6

disp(' ------ interpolation order-6-polynomial -----')
order = 6;
% coefficients of sampled coordinates
Xmatrix = [ones(M,1),Mset,Mset.^2,Mset.^3,Mset.^4,Mset.^5,Mset.^6];
% estimate params by regression
estpar  = regress(ys,Xmatrix(:,1:order+1));

% reconstruction
xc=(1:N)';
Xpred=[ones(N,1),xc,xc.^2,xc.^3,xc.^4,xc.^5,xc.^6];
xpred = Xpred(:,1:order+1)*estpar;

subplot(4,3,6); hold on;
plot(1:N,xpred,'-r','LineWidth',2)
plot(Mset,ys,'.b','MarkerSize',15)
plot(1:N,x,'--k','LineWidth',1)
ylim([ymin,ymax])
title(strcat('best order-',num2str(order),'-polynomial'))

%% interpolation with very low smoothness

disp(' ------ interpolation very high smoothness -----')
sigma_e = 0.02;
xest = estimate_profile_robust...
    (N,select,ys,sigma_e,1,1,0,[0,0,0,0],0,0);

subplot(4,3,7);hold on;
plot(1:N,xest,'-r','LineWidth',2)
plot(Mset,ys,'.b','MarkerSize',15)
plot(1:N,x,'--k','LineWidth',1)
ylim([ymin,ymax])
title(['smooth reconstruction, $\sigma_e =',num2str(sigma_e),'$'])

%% interpolation with medium smoothness

disp(' ------ interpolation medium smoothness  -----')
sigma_e = 0.1;
xest = estimate_profile_robust...
    (N,select,ys,sigma_e,1,1,0,[0,0,0,0],0,0);

subplot(4,3,8);hold on;
plot(1:N,xest,'-r','LineWidth',2)
plot(Mset,ys,'.b','MarkerSize',15)
plot(1:N,x,'--k','LineWidth',1)
ylim([ymin,ymax])
title(['smooth reconstruction, $\sigma_e =',num2str(sigma_e),'$'])

%% interpolation with very high smoothness

disp(' ------ interpolation high smoothness -----')
sigma_e = 0.5;
xest = estimate_profile_robust...
    (N,select,ys,sigma_e,1,1,0,[0,0,0,0],0,0);

subplot(4,3,9);hold on;
plot(1:N,xest,'-r','LineWidth',2)
plot(Mset,ys,'.b','MarkerSize',15)
plot(1:N,x,'--k','LineWidth',1)
ylim([ymin,ymax])
title(['smooth reconstruction, $\sigma_e =',num2str(sigma_e),'$'])


%% interpolation medium flatness

disp(' ------ interpolation high flatness -----')
sigma_e = 0.2;
xest = estimate_profile_robust_flat...
    (N,select,ys,sigma_e,1,1,0,[0,0,0,0],0);

subplot(4,3,10);hold on;
plot(1:N,xest,'-r','LineWidth',2)
plot(Mset,ys,'.b','MarkerSize',15)
plot(1:N,x,'--k','LineWidth',1)
ylim([ymin,ymax])
title(['flat reconstruction, $\sigma_e =',num2str(sigma_e),'$'])


%% interpolation with high flatness

disp(' ------ interpolation low flatness -----')
sigma_e = 1.0;
xest = estimate_profile_robust_flat...
    (N,select,ys,sigma_e,1,1,0,[0,0,0,0],0);

subplot(4,3,11);hold on;
plot(1:N,xest,'-r','LineWidth',2)
plot(Mset,ys,'.b','MarkerSize',15)
plot(1:N,x,'--k','LineWidth',1)
ylim([ymin,ymax])
title(['flat reconstruction, $\sigma_e =',num2str(sigma_e),'$'])

%% interpolation mixed

disp(' ------ interpolation left: flat / right: smooth -----')
sigma_e1 = 0.2;
sigma_e2 = 0.5;
fs = [ones(N1,1);ones(N-N1,1)*2];
xest = estimate_profile_robust_flat_smooth...
    (N,fs,select,ys,sigma_e1,sigma_e2,1,1,0,[0,0,0,0],0);

subplot(4,3,12);hold on;
plot(1:N,xest,'-r','LineWidth',2)
plot(Mset,ys,'.b','MarkerSize',15)
plot(1:N,x,'--k','LineWidth',1)
ylim([ymin,ymax])
title(strcat('mixed reconstruction'))

%% three samples
disp(' --------------------------------------------------')
disp(' --------- Fig. 16.4 Three sample profiles --------')
disp(' --------------------------------------------------')

N = 200;
N1 = 70;
N2 = N-N1;

figure('name','Fig. 16.4: Three sample profiles','color','w',...
    'Position',[0.3*ss(1),0.2*ss(2),0.5*ss(1),0.3*ss(2)]);

for samples = 1:3
    
    [x1,y1] = generate_observed_ARp(N1,1,1,sigma_e,0.5);
    % enforce slope 0
    for i=1:N1
        y1(i)=y1(i) -  (i-1)*(x1(N1)-x1(1))/(N1-1);
        x1(i)=x1(i) -  (i-1)*(x1(N1)-x1(1))/(N1-1);
    end
    
    [x2,y2] = generate_observed_ARp(N2,2,1,sigma_n,0.5);
    % enforce slope 0
    for i=1:N2
        y2(i)=y2(i) -  (i-1)*(x2(N2)-x2(1))/(N2-1);
        x2(i)=x2(i) -  (i-1)*(x2(N2)-x2(1))/(N2-1);
    end
    
    x = [x1;x2+x1(N1)];
    y = [y1;y2+y1(N1)];

    ymin = -250;
    ymax = 250;

    subplot(1,3,samples)
    hold on
    plot(1:N,x,'-k','LineWidth',2)
    ylim([ymin,ymax])
    title(['sample profile no: ',num2str(samples)])
end
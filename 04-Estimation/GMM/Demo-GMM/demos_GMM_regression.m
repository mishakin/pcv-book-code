%% test GMM with linear regression with check
%
% 1. Single example
% 2. checking the correctness
% 3. the effect of varying noise level
% 4. the effect of varying density of observations
% 
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 

close all
clear sigma_x12

addpath(genpath('../../../General-Functions/'));
addpath(genpath('../Functions-GMM/'));


disp('============================================')
disp('============================================')
disp('             Linear regression              ')
disp('            l = x[1] + x[2] * t             ')
disp('--------------------------------------------')

%% General control parameters

% choose one or several of the following four parts of the example
%choice=[1,0,0,0];  % single example
%choice=[0,1,0,0];  % checking the correctness%
%choice=[0,0,1,0];  % the effect of varying noise level
%choice=[0,0,0,1];  % the effect of varying density of observations
choice=[1,1,1,1];

% initialization of random numbers
% = 0 CPU-time dependent
% > 0 fixed
%init_rand  = 0;
init_rand   = 15;     % standard = 15

%% random number intialization
init_rand_seed(init_rand);

%%    
% unknown parameters
xt  = [0.5,1]';
U   = 2;
disp(['True_parameters = [', num2str(xt'),']'])

% Significance level
S = 0.99;

% non-centrality parameter 
delta0 = 4;

%% Individual control parameters

% 1. for diagnostics: set of parameters to be evaluated
%------------------------------------------------------
rU = 1;  % can be [], 1, 2, or [1,2]; 1 for offset, 2 for scale


% 2. for checking the result
% --------------------------
% type of generated data 
% data_type = -10;  % 10 random position in 1:10
data_type2 =   1;  % 4 observations with leverage point at i=4;

% number of samples for estimation
K2 = 25;


% 3. for varying sigma 
% --------------------
% type of generated data
% data_type = -10;  % 10 random position in 1:10
data_type3 =   1;  % 4 observations with leverage point at i=4;

% varying sigmas
sigmas = [0.01,0.05,0.1,0.15,0.20];

% number of samples per sigma
K3 = 25;


% 4. for varying the density/number N of observations
% ---------------------------------------------------

% varying N's
%Ns = [5,10,20,40,80,160,320];
Ns = round(exp(1+2*log(2):0.5*log(2):7)/exp(1));

% number of samples per N
K4 = 5;




%% 1. Single example #########################################################
if choice(1) == 1
    disp('============================================')
    disp('         Example with evaluation            ')
    disp('............................................')

    % initializing random numbers-------------------------------------------
    
    init_rand_seed(init_rand);
    
    
    if ~isempty(rU)
        if rU(1) == 1
            disp('parameters_to_be_evaluated = intercept')
        else
            disp('parameters_to_be_evaluated = slope')
        end
    end
    
    % times of observation ------------------------------------------------
    tv    = [-1,1,2,14]';
    N   = length(tv);
    disp(['Number of observations = ', num2str(N)])%
    Am  = [ones(N,1), tv];
    av  = zeros(N,1);
    
    % stochastical model --------------------------------------------------
    sigma= 0.5;
    Cov_ll = sigma^2 * eye(N);
    
    % observations --------------------------------------------------------
    true_error =randn(N,1)*sigma;
    lv = Am*xt+true_error;
    % add outlier err(4)=+4 in 4th observation
    lv(4)=lv(4)-4;
    
    disp('    no         t        tl        te         l')
    disp([(1:4)',tv,Am*xt,lv-Am*xt,lv ])
   
    % estimate parameters -------------------------------------------------
    start=cputime;
    [est_x,Cov_xx,sigma_0q,R,vv,zv,riv,nabla_lv,muv,muv1,uv1q,uv2]=...
        GaussMarkovModelLinear(lv,Cov_ll,Am,av,rU);
    CPU_time = (cputime-start)                                              %#ok<*NOPTS>
    
    warning('off','all');
    
    % --------------- output and diagnostics ------------------------------
    ScrS = plot_init;  
    f1 = figure('Name','Regression and Diagnostic','Color','w','Position',[20 ScrS(2)-0.7*ScrS(2)-100 0.6*ScrS(1) 0.7*ScrS(2)]);
    subplot(2,4,1)
    hold on
    plot(tv,lv,'.r','Markersize',12)       % given data
    plot([min(tv),max(tv)],[[1,min(tv)]*est_x,[1,max(tv)]*est_x],'-b')
    xlim([min(tv)-1,max(tv)+1])
    title({'linear regression', strcat('$\sigma_{l} = ',num2str(sigma),'$')})
    
    xlabel('t')
    
    
    subplot(2,4,2)
    hold on
    plot(tv,vv,'.r')
    plot([min(tv),max(tv)],[0,0],'-b')
    title({'residuals', strcat('max($v_i$) = ',num2str(max(abs(vv))))})
    xlim([min(tv)-1,max(tv)+1])
    
    xlabel('t')
    
    subplot(2,4,3)
    hold on
    plot(tv,riv,'.r')
    title({'redundancy number', strcat('min($r_i$) = ',num2str(min(riv)))})
    xlim([min(tv)-1,max(tv)+1])
    
    xlabel('t')
    
    subplot(2,4,4)
    hold on
    est_e=-vv./riv;
    plot(tv,est_e,'.r')
    plot([min(tv),max(tv)],[0,0],'-b')
    title({'estimated error', strcat('$\max(|\hat{e}_i|)) = $',num2str(max(abs(est_e))))})
    xlim([min(tv)-1,max(tv)+1])
    ylim([min(est_e)-1,max(est_e)+1])
    xlabel('t')
    
    
    
    subplot(2,4,5)
    hold on
    plot(tv,nabla_lv,'.r')
    title({'min detectable outliers', strcat('$\max(\nabla_{{0l}_i}) = $',...
        num2str(max(nabla_lv)))})
    xlim([min(tv)-1,max(tv)+1])
    xlabel('t')
    
    
    subplot(2,4,6)
    hold on
    plot(tv,zv,'.r')
    plot([min(tv),max(tv)],[2.58,2.58],'-r')
    plot([min(tv),max(tv)],-[2.58,2.58],'-r')
    title({'test statistics', strcat('max($z_i$) = ',num2str(max(abs(zv))))})
    xlim([min(tv)-1,max(tv)+1])
    xlabel('t')
    
    
    
    subplot(2,4,7)
    hold on
    plot(tv,muv,'.r')
    title({'sensitivity factors', strcat('max($\mu_i$) = ',num2str(max(muv)))})
    xlim([min(tv)-1,max(tv)+1])
    xlabel('t')
    
    
    %
    
    disp('........................................................')
    disp('                      diagnostics                       ')
    disp('........................................................')
    disp(strcat('    no         t         error     v         r     est_error',...
            '     z         z*   nabla_0l  nabla_0^*l' ));
    o=[(1:N)', tv, lv-Am*xt, vv, riv, -vv./riv, ...
        zv, zv.*sqrt(riv), delta0*sigma./sqrt(riv), delta0*sigma./riv ]
    
    Estimated_parameters      = est_x'
    Theoretical_precision     = sqrt(diag(Cov_xx))'
    Redundancy                = R
    Estimated_sigma_0         = sqrt(sigma_0q)
    
    [m,i]=max(abs(true_error));
    disp(horzcat('Maximal true error                 = ',num2str(m,'% 12.5f'),...
        ' at: ', num2str(i,'% 3.0f')));
    [m,i]=max(abs(vv));
    disp(horzcat('Maximal residual                   = ',num2str(m,'% 12.5f'),...
        ' at: ', num2str(i,'% 3.0f')));
    [m,i]=min(riv);
    disp(horzcat('Minimal redundancy number          = ',num2str(m,'% 12.5f'),...
        ' at: ', num2str(i,'% 3.0f')));
    
    [m,i]=max(abs(zv));
    if max(abs(zv)) > norminv(S,0,1)
        disp(horzcat('Maximal test statistic             = ',num2str(m,'% 12.5f'),...
            ' at: ', num2str(i,'% 3.0f')), ' ***');
    else
        disp(horzcat('Maximal test statistic             = ',num2str(m,'% 12.5f'),...
            ' at: ', num2str(i,'% 3.0f')));
    end
    
    [m,i]=max(nabla_lv);
    disp(horzcat('Max. of minimal detectable outlier = ',num2str(m,'% 12.5f'),...
        ' at: ', num2str(i,'% 3.0f')));
    
    [m,i]=max(muv);
    disp(horzcat('Max. sensitivity factor            = ',num2str(m,'% 12.5f'),...
        ' at: ', num2str(i,'% 3.0f')));
    
    if ~isempty(rU)
        [m,i]=max(muv1);
        disp(horzcat('Max. sensitivity factor selected   = ',num2str(m,'% 12.5f'),...
            ' at: ', num2str(i,'% 3.0f')));
    end
    
    
    % show sensitivity factors --------------------------------------------
    m=2.995;range=0:m/200:m;
    N=10;
    mux=sqrt((1+range.^2)./(N-(1+range.^2)));
    mu1=sqrt(1./(N-(1+range.^2)));
    mu2=sqrt(range.^2./(N-(1+range.^2)));
    
    f2 = figure('Name','Sensitivity factor mu_x','Color','w','Position',[30+0.35*ScrS(1) ScrS(2)-0.4*ScrS(2)-110 0.3*ScrS(1) 0.4*ScrS(2)]);
    plot(range,mux,'-r','LineWidth',2);
    title('sensitivity factor $\mu_x$ = f(distance)')
    xlabel('$d$')
    ylabel('$\mu_x$')
    ylim([0,5])
    
    f3 = figure('Name','Sensitivity factor mu_offset','Color','w','Position',[0.7*ScrS(1) ScrS(2)-0.4*ScrS(2)-110 0.3*ScrS(1) 0.4*ScrS(2)]);
    hold on
    plot(range,mu1,'-r','LineWidth',2);
    plot(range,mux,'-g','LineWidth',2);
    title('sensitivity factor $\mu_{\mbox{offset}}$ (red) and $\mu_x$ (green) = f(distance)')
    ylim([0,5])
    xlabel('$d$')
    ylabel('$\mu_1,\ \mu_x$')
    ylim([0,5])
    
    f4 = figure('Name','Sensitivity factor mu_offset','Color','w','Position',[30 20 0.3*ScrS(1) 0.4*ScrS(2)]);
    hold on
    plot(range,mu2,'-r','LineWidth',2);
    plot(range,mux,'-g','LineWidth',2);
    title('sensitivity factor $\mu_{scale}$ (red) and $\mu_x$ (green) = f(distance)')
    ylim([0,5])
    xlabel('$d$')
    ylabel('$\mu_2,\ \mu_x$       ')
    ylim([0,5])
    
end



%% 2. checking the correctness
if choice(2) == 1
    
    disp('========================================================')
    disp('             checking the correctness                   ')
    disp('........................................................')
    
    % initializing random numbers -----------------------------------------
    % the next two lines can be exchanged
    init_rand_seed(init_rand)
      
    % K: number of samples ------------------------------------------------
    K=K2;  % relative accuracy is sqrt(1/K)
    
    
    switch data_type2
        case 1
            tv = [-1,1,2,14]';
        otherwise
            tv = rand(10,1)*10;
            tv = sort(tv);
            % sigma  = 1.5;
    end
    
    N   = length(tv);
    disp(['Number of observations = ', num2str(N)])
    Am  = [ones(N,1), tv];
    av  = zeros(N,1);
    
    % stochastical model --------------------------------------------------
    sigma= 0.5;
    Cov_ll = sigma^2 * eye(N);
    
    % K samples
    est_xs   = zeros(K,2);
    est_s0qs = zeros(K,1);
    start=cputime;
    for k = 1:K
        true_error = randn(N,1)*sigma;
        lv = Am*xt+true_error;
        
        % estimate parameters
        
        [est_x,Cov_xx,sigma_0q,R,vv,zv,riv,nabla_lv,muv,muv1,uv1q,uv2]=...
            GaussMarkovModelLinear(lv,Cov_ll,Am,av,rU);
        est_s0qs(k)   = sigma_0q;
        est_xs(k,:)   = est_x';
    end
    
    CPU_time_per_sample=(cputime-start)/K
    
    warning('off','all');
    
    % output and diagnostics ----------------------------------------------
       
    f5 = check_estimation_result(R,xt,Cov_xx,est_s0qs,est_xs,S,'linear regression');
    set(f5,'Name','Variance Factors','Position',[30+0.35*ScrS(1) 20 0.3*ScrS(1) 0.4*ScrS(2)])
    
end

%% 3. the effect of varying noise level

if choice(3) == 1
    clear tv 
    disp('========================================================') 
    disp('          effect of varying noise level                 ')
    disp('........................................................')
    % initializing random numbers -----------------------------------------
    init_rand_seed(init_rand)
    
    switch data_type3
        case 1
            tv    = [-1,1,2,14]';
        otherwise
            tv=rand(10,1)*10;
            tv=sort(tv);
            sigma  = 1.5;
    end
    
    N   = length(tv);
    disp(['Number of observations = ', num2str(N)])
    Am  = [ones(N,1), tv];
    av  = zeros(N,1);
    
    
    K = K3;
    disp(['Number of samples per sigma = ', num2str(K)])
    
    
    disp(['sigma''s = [', num2str(sigmas),']'])
    L = length(sigmas);
    
    % for all sigmas
    sigma_x12 = zeros(L,2);
    for l=1:L
        sigma= sigmas(l);
        Cov_ll = sigma^2 * eye(N);
        %% generate K samples and estimate
        est_xs   = zeros(K,2);
        est_s0qs = zeros(K,1);
        
        % for all samples
        start=cputime;
        for k=1:K
            
            true_error = randn(N,1)*sigma;
            lv = Am*xt+true_error;
            
            % estimate parameters
            
            [est_x,Cov_xx,sigma_0q,R,vv,zv,riv,nabla_lv,muv,muv1,uv1q,uv2]=...
                GaussMarkovModelLinear(lv,Cov_ll,Am,av,rU);
            est_s0qs(k)   = sigma_0q;
            est_xs(k,:)   = est_x';
        end
        
        warning('off','all');
        
        Cov_est_x_emp = cov(est_xs);
        %  set sigmas
        sigma_x12(l,:)=...
            [sqrt(Cov_est_x_emp(1,1)),sqrt(Cov_est_x_emp(2,2))];
        
    end
    
    % plot dependencies
    f6 = figure('Name','Effect of varying noise level ','Color','w','Position',[30 20 0.6*ScrS(1) 0.7*ScrS(2)]);
    subplot(2,2,1)
    hold on
    plot(sigmas,sigmas/sigmas(L)*sqrt(Cov_xx(1,1))','--b','Linewidth',2);
    plot(sigmas,sigma_x12(:,1)','-r','Linewidth',2);
    plot(sigmas,sigma_x12(:,1)','.k','MarkerSize',15);
    title('$\sigma_1(\sigma)$')
    xlabel('$\sigma$')
    subplot(2,2,2)
    hold on
    plot(sigmas,sigmas/sigmas(L)*sqrt(Cov_xx(2,2))','--b','Linewidth',2);
    plot(sigmas,sigma_x12(:,2)','-r','Linewidth',2);
    plot(sigmas,sigma_x12(:,2)','.k','MarkerSize',15);
    title('$\sigma_2(\sigma)$')
    xlabel('$\sigma$')
    subplot(2,2,3)
    hold on
    plot(sigmas,ones(1,L)./sigmas(L)*sqrt(Cov_xx(1,1))','--b','Linewidth',2);
    plot(sigmas,sigma_x12(:,1)./sigmas','-r','Linewidth',2);
    plot(sigmas,sigma_x12(:,1)./sigmas','.k','MarkerSize',15);
    title('$\sigma_1(\sigma)/\sigma$')
    ylim([0,1.2*max(sigma_x12(:,1)./sigmas')])
    xlabel('$\sigma$')
    subplot(2,2,4)
    hold on   
    plot(sigmas,ones(1,L)./sigmas(L)*sqrt(Cov_xx(2,2))','--b','Linewidth',2);
    plot(sigmas,sigma_x12(:,2)./sigmas','-r','Linewidth',2);
    plot(sigmas,sigma_x12(:,2)./sigmas','.k','MarkerSize',15);
    title('$\sigma_2(\sigma)/\sigma$')
    ylim([0,1.2*max(sigma_x12(:,2)./sigmas')])
    xlabel('$\sigma$')
    
end

%% 4. the effect of varying density of observations

if choice(4) ==1
    
    disp('========================================================')
    disp('     effect of varying density of observations          ')
    disp('........................................................')
    
    % initializing random numbers -----------------------------------------
    % the next two lines can be exchanged
    init_rand_seed(init_rand)
    
    disp(strcat('Ns = [',num2str(Ns),']'))

    L = length(Ns);
    K = K4;
    
    disp(['Number of samples per N and estimates = ', num2str(K)])

    % sigma
    sigma=0.02;
    
    % needs to be set ...
    rU=1;
    
    lm=0;
    % for all Ns
    for l=1:L
        N=Ns(l);
        
        % also repeat sqrt(K) times Ns(l)
        v1=0;
        v2=0;
        for m=1:K
            
            % sorted observations in interval [0..100]
            %tv=(rand(N,1))*100;
            tv=(1/(2*N):100/(N):100-1/(2*N))';
            tv = sort(tv);                                                  %#ok<TRSRT>
            Am  = [ones(N,1), tv];
            av  = zeros(N,1);
            
            Cov_ll = sigma^2 * eye(N);
            
            % generate K samples and estimate
            est_xs   = zeros(K,2);
            est_s0qs = zeros(K,1);
            
            % for all samples
            start=cputime;
            for k=1:K
                
                true_error = randn(N,1)*sigma;
                lv         = Am*xt + true_error;
                
                % estimate parameters
                
                [est_x,Cov_xx,sigma_0q,R,vv,zv,riv,nabla_lv,muv,muv1,uv1q,uv2]=...
                    GaussMarkovModelLinear(lv,Cov_ll,Am,av,rU);
                est_s0qs(k)   = sigma_0q;
                est_xs(k,:)   = est_x';
            end
            
            warning('off','all');
            
            Cov_est_x_emp = cov(est_xs);
            %  set sigmas
            v1=v1+Cov_est_x_emp(1,1);
            v2=v2+Cov_est_x_emp(2,2);
            
        end
        sigma_x12(l,:)=sqrt([v1/K,v2/K]);
    end
    
    % plot dependencies
    f7 = figure('Name','Effect of varying density of observations ','Color','w','Position',[ScrS(1)-0.6*ScrS(1) 20 0.6*ScrS(1) 0.7*ScrS(2)]);
    subplot(2,2,1)
    hold on
    plot(Ns,sqrt(Cov_xx(1,1))./sqrt(Ns/Ns(L))','--b','Linewidth',2);
    plot(Ns,sigma_x12(:,1)','-r','Linewidth',2);
    plot(Ns,sigma_x12(:,1)','.k','MarkerSize',15);
    title('$\sigma_1(N)$')
    xlabel('$N$')
    subplot(2,2,2)
    hold on
    plot(Ns,sqrt(Cov_xx(2,2))./sqrt(Ns/Ns(L))','--b','Linewidth',2);
    plot(Ns,sigma_x12(:,2)','-r','Linewidth',2);
    plot(Ns,sigma_x12(:,2)','.k','MarkerSize',15);
    title('$\sigma_2(N)$')
    xlabel('N')
    subplot(2,2,3)
    hold on
    plot(Ns,sqrt(Cov_xx(1,1))./sqrt(1/Ns(L))*ones(1,L),'--b','Linewidth',2);
    plot(Ns,sqrt(Ns').*sigma_x12(:,1),'-r','Linewidth',2);
    plot(Ns,sqrt(Ns').*sigma_x12(:,1),'.k','MarkerSize',15);
    title('$\sqrt{N} \sigma_1(N)$')
    xlabel('$N$')
    ylim([0,1.2*max(sqrt(Ns').*sigma_x12(:,1))])
    subplot(2,2,4)
    xlabel('$N$')
    hold on
    plot(Ns,sqrt(Cov_xx(2,2))./sqrt(1/Ns(L))*ones(1,L),'--b','Linewidth',2);
    plot(Ns,sqrt(Ns').*sigma_x12(:,2),'-r','Linewidth',2);
    plot(Ns,sqrt(Ns').*sigma_x12(:,2),'.k','MarkerSize',15);
    title('$\sqrt{N} \sigma_2(N)$')
    ylim([0,1.2*max(sqrt(Ns').*sigma_x12(:,2))])
    xlabel('$N$')
    
end

figure(f1)
figure(f2)
figure(f3)
figure(f6)
figure(f7)
figure(f4)
figure(f5)
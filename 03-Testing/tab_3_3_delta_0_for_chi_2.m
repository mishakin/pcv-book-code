%% Table 3.3 to 3.5 page 68
% determine delta_0 for chi^2 distribution
%
% Wolfgang Förstner 2015
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de


disp('---- determine lower bounds delta for non-central chi-distribution ----')

%% set parameters

% fix alpha's
alphas = [0.003:-0.00005:0.002];
Nalpha = length(alphas);

% set beta's
betas = [0.8] %,0.9,0.95];

% fix degrees of freedom
dof =[2]; %,3,4,8,10,20,50];
Ndof = length(dof);


%% complete table including left and top margin
table_alpha_d = zeros(1+Nalpha,1+Ndof);
table_alpha_d(2:Nalpha+1,1) = alphas;
table_alpha_d(1,2:1+Ndof) = dof;

% for all beta's
for ib = 1:length(betas);

    beta_0 = betas(ib);
    
    % for all alpha's
    for ia = 1:Nalpha
        
        alpha = alphas(ia);
        
        % for all dof's
        for idf = 1:length(dof)
        
            df = dof(idf);
            
            % critical value for chi-square test
            c = chi2inv(1-alpha,df);
            
            % set boundaries for searching for delta_0
            mindelta = 0.1;
            maxdelta = 6*df;
            meandelta = (maxdelta+mindelta)/2;
            
            % use noncentral chi-square cumulative distribution function to
            % determine boundaries
            minpower = 1 - ncx2cdf(c,df,mindelta^2);
            maxpower = 1 - ncx2cdf(c,df,maxdelta^2);
            
            % search binary, half intervals until right and left interval
            % ends have the same power
            while maxpower-minpower > 10^-10
                
                d_power=maxpower-minpower;
                meanpower= 1-ncx2cdf(c,df,meandelta^2);
                
                if meanpower > beta_0
                    maxdelta  =  meandelta;
                    maxpower = 1-ncx2cdf(c,df,maxdelta^2);
                    meandelta = (maxdelta+mindelta)/2;
                else
                    mindelta = meandelta;
                    minpower = 1-ncx2cdf(c,df,mindelta^2);
                    meandelta = (maxdelta+mindelta)/2;
                end                
                
            end
            
            table_alpha_d(ia+1,idf+1) = meandelta;
            
        end
    end
    beta = beta_0 %#ok<*NOPTS>
    t_alpha_0_d = table_alpha_d

end
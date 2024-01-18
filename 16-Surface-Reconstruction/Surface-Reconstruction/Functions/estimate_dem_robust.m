%% robustly estimate correction to dem
%
% Optimization function using curvature for regularization
%   energy = sum_m rho(l_m-a_m)/sigma_n) + 
%             1/sigma_d sum_ij d_ii^2+2d_ij^2+djj^2
% l_m = observed heights, interpolated bilinearly
% a_m unknown grid heights, d_ij second derivatives
%
% [A,weights,dx,v,W,Nmatrix]=...
%   estimate_dem_robust...
%    (N,U,Np,Nr,Mc,poi,xa,weights,sigma_k,iter,iter_switch,type_robust,L1,...
%    out_in,print_type,plot_type)
%
% N           = number of observations (points+second derivatives)
% U           = number of unknowns (dem)
% Np          = number of points
% Nr          = number of rows
% Mc          = number oc columns 
% poi         = Npx4 matrix of points (x,y,h,sigma), x,y, referring to grid
% xa          = approximate values for dem
% weights     = Nx1 vector of weights in [0,1]
% sigma_k     = standard deviation of second derivatives
% iter        = current iteration for controling weighting
% iter_switch = iteration number where to switch type of weighting
% type_robust = 0,1,2,3, (00,01,10,11) for points and dem
% L1          = 0/1 choice of weight function (0 = exponential, 1 = L1)
% out_in      = boolean Np-vector indicating inliers
% print_type  = 0,1,2
% plot_type   = 0,1,2
%
% A         = NxU sparse design matrix
% weights   = new weights
% dx        = delta dem
% v         = normalized residuals
% W         = weigths of obs. combining uncertainty and robust factor
% Nmatrix   = normal equation matrix
%
% wf 7/2014, 2/2018
%

function [A,weights,dx,v,W,Nmatrix]=...
    estimate_dem_robust...
    (N,U,Np,Nr,Mc,poi,xa,weights,sigma_k,iter,iter_switch,type_robust,L1,...
    out_in,print_type,plot_type);
    
typec=1;

%% factor for update weights indicating outliers
k= 2;

%% start estimation
% initiate size of A
            if print_type > 0
                start_time_sol = cputime;
            end
A = spalloc(N,U,6*N);
b = zeros(N,1); 
W = spalloc(N,N,N);
Ac = spalloc(N,U,6*N);
bc = zeros(N,1); 
Wc = spalloc(N,N,N);

        % observations, normalized with 1/sqrt of variances
        tic
        for p = 1:Np
            x = poi(p,1);
            y = poi(p,2);
            w = 1/poi(p,4)*sqrt(weights(p))*out_in(p);
            i = floor(x);
            j = floor(y);
            aa = x-i;
            bb = y-j;
            maa = 1-aa;
            mbb=  1-bb;
            g1=maa*mbb;
            g2=aa*mbb;
            g3=maa*bb;
            g4=aa*bb;
%             %delta_x
%             A(p,i+j*Nr+1)       = (1-aa)*(1-bb)*w;
%             A(p,i+j*Nr+2)       =    aa *(1-bb)*w;
%             
%             A(p,i+(j+1)*Nr+1)   = (1-aa)* bb   *w;
%             A(p,i+(j+1)*Nr+2)   =    aa * bb   *w;
%             delta_l =  poi(p,3) - ...
%                 ((1-aa)*(1-bb)*xa(i+j*Nr+1)    +...
%                 aa *(1-bb)*xa(i+j*Nr+2)    +...
%                 (1-aa)* bb   *xa(i+(j+1)*Nr+1)+...
%                 aa * bb   *xa(i+(j+1)*Nr+2));
            A(p,i+j*Nr+1)       = g1*w;
            A(p,i+j*Nr+2)       = g2*w;
            A(p,i+(j+1)*Nr+1)   = g3*w;
            A(p,i+(j+1)*Nr+2)   = g4*w;
            delta_l =  poi(p,3) - ...
                (g1*xa(i+j*Nr+1)    +...
                 g2*xa(i+j*Nr+2)    +...
                 g3*xa(i+(j+1)*Nr+1)+...
                 g4*xa(i+(j+1)*Nr+2));
            b(p)                =   delta_l * w;
            W(p,p)=w;
        end
        
        time_for_A_points=toc


            if print_type > 0
                time_for_building_An=toc
            end

%% priors
tic
Icc = Np;                    % first index for column curvatures
Irr = Np + (Nr-2)*Mc;        % first index for row curvature
Irc = Np + (Nr-2)*Mc + Nr*(Mc-2);    % first index for torsion 
% coefficients for column curvature
        % cc  1 
        %    -2 
        %     1 
Jcc = [-1,0,+1];  % index vector
Vcc = [ 1 -2 1];                    % coefficients
wcc = 1/sigma_k;%/sqrt(6);  % /norm(Vcc)
% coefficients for row curvature
        % rr  1 -2  1
Jrr = [-Nr,0,+Nr];  % index vector
Vrr = [1 -2  1];                    % coefficients
wrr = 1/sigma_k;% /sqrt(6); % /norm(Vrr)
% coefficients for torsion
        % rc  1  -1 
        %    -1   1
%                  ^
Jrc = [-Nr-1,-Nr,-1,0];      % indices
Vrc = [1 -1 -1 1];                % coefficients
%wrc = 1/sigma_k*2;  %/sqrt(2); %/2; % /norm(Vrc)
wrc = 1/sigma_k/sqrt(2);      % quadratic variation
% built up of priors in A
% column curvature
for j=1:Mc
    for i=2:Nr-1
        IU  = i+(j-1)*Nr;
        %[i,j,IU]
        Icc = Icc+1;
        A(Icc,IU+Jcc) = Vcc * wcc * sqrt(weights(Icc));
        delta_cc      = 0-xa(IU+Jcc)' * Vcc';
        b(Icc)        = delta_cc * wcc * sqrt(weights(Icc));
        W(Icc,Icc)    = wcc * sqrt(weights(Icc));
    end
end

% row curvature
for j=2:Mc-1
    for i=1:Nr
        IU  = i+(j-1)*Nr;
        %[i,j,IU]
        Irr = Irr+1;
        A(Irr,IU+Jrr) = Vrr * wrr * sqrt(weights(Irr));
        delta_rr      = 0-xa(IU+Jrr)' * Vrr';
        b(Irr)        = delta_rr * wrr * sqrt(weights(Irr));
        W(Irr,Irr)    = wrr * sqrt(weights(Irr));
    end
end



% torsion
for j=2:Mc
    for i=2:Nr
        IU  = i+(j-1)*Nr;
        %[i,j,IU]
        Irc = Irc+1;
        A(Irc,IU+Jrc) = Vrc * wrc * sqrt(weights(Irc));
        delta_rc      = 0-xa(IU+Jrc)' * Vrc';
        b(Irc)        = delta_rc * wrc * sqrt(weights(Irc));
        W(Irc,Irc)    = wrc * sqrt(weights(Irc));
    end
end
time_for_A_prior=toc
%
            if print_type > 1
                if U < 1000
                    Ab = 2*[full(A)]
                    null_A_prior = U-rank(full(A(U+1:end,:)));
                    rankA = rank(full(A(U+1:end,:)));
                end


            end
%time_for_building_Ae=toc
%% solve

% sort
tic
Nmatrix  = A'*A;


            if iter==1 && plot_type > 1
                size_N_matrix = size(Nmatrix);
                non_zeros_N_matrix =nnz(Nmatrix);
                figure
                subplot(1,3,1)
                spy(Nmatrix)
                title('N')
            end

permx    = symamd(Nmatrix);
Aperm    = A(:,permx);




            if iter==1 && plot_type > 1
                subplot(1,3,2)
                spy(Aperm'*Aperm)
                title('N, sorted')
            end


[C,R]      = qr(Aperm,b,0);

dxperm   = R\C;
Pm       = speye(U);
Pm       = Pm(:,permx);
dx       = Pm*dxperm;



                time_for_QR_sorted=toc

            if iter==1 && plot_type > 1
                subplot(1,3,3)
                spy(R)
                fill_in = nnz(R)-(nnz(Nmatrix)/2+U/2)
                percentage_nnz   = nnz(R)/(U*(U+1)/2)
                relative_fill_in = fill_in/(nnz(Nmatrix)/2+U/2)
                title(strcat('reduced N sorted, fill in =',num2str(fill_in)))
            end


            if iter==1 && plot_type > 1
                figure
                subplot(1,2,1)
                spy(A)
                title('A')
                subplot(1,2,2)
                spy(Aperm)
                title('A sorted')
            end

dx_sorted = dx;

            if U < 1600 && plot_type > 0
                Aqt = R\Aperm';
                figure
                spy(Aqt')
                title('R\A sorted')
            end

            if print_type > 0
                time_for_solution=cputime-start_time_sol
            end

% residuals (normalized)
v = A*dx-b;



% estimated variance factor
eso = norm(v)/sqrt(N-U);

if type_robust > 0
    
    
    % original non-normalized residuals
    vori  = v;
    vori(1:Np)     = (v(1:Np).*poi(1:Np,4)./sqrt(weights(1:Np)));
    vori(Np+1,end) = (v(Np+1,end)*sigma_k./sqrt(weights(Np+1,end)));
    vs=vori;
    %vs=v;
        
    % sigma_0 for points
    eso_p = median(abs(v(1:Np)))*1.48;
    %es_p  = median(abs(vs(1:Np)))*1.48;
    
    % sigma_0 for curvatures
    eso_k = median(abs(v(Np+1:end)))*1.48;
    %es_k  = median(abs(vs(Np+1:end)))*1.48;
                if print_type > 0
                    est_s0_s0p_s0k=[eso,eso_p,eso_k]
                end

    check_AtV_null=norm((A'*v));
    max_v_points = max(abs(v(1:Np)));
    max_v_curvat = max(abs(vs(Np+1,end)));

    % critical values
    kp = k*eso_p;
    kk = k*eso_k;
%     kp = k*es_p;
%     kk = k*es_k;
   
    %kp = k;
    
    
                if print_type > 0
                    criticalvalue_kp_sigmap_kk_sigma_k=[kp,poi(1,4),kk,sigma_k]
                end
end

if iter == 1 && type_robust > 0 
    weights=ones(N,1);
else
    if type_robust > 0 && L1
        % dem
        if type_robust == 1 || (type_robust == 3 && iter > 2*iter_switch)
            %weights(Np+1:end) = min(1,1./(abs(vs(Np+1:end))/kk+eps));
            weights(Np+1:end) = min(1,1./(abs(v(Np+1:end))/kk+eps));
        end
        % points
        if type_robust == 2 || (type_robust ==3 && iter <= iter_switch)
            weights(1:Np) = min(1,1./(abs(v(1:Np))/kp+eps)).*out_in(1:Np);
        end
    else
        if type_robust == 1 || (type_robust == 3 && iter > 2*iter_switch)
            %weights(Np+1:end) = exp(-abs(vs(Np+1:end)).^2/kk^2/2)+0.0001;
            weights(Np+1:end) = exp(-abs(v(Np+1:end)).^2/kk^2/2)+0.0001;
            %weights(Np+1:end) = exp(-abs(vs(Np+1:end))./kk+0.0001);
        end
        if type_robust == 2 || (type_robust ==3 && iter <= iter_switch)          
            weights(1:Np) = (exp(-abs(v(1:Np)).^2/kp^2/2)+0.0001).*out_in(1:Np);
        end
    end
end
%iter=iter
abs_dx = norm(dx);
min_w_p=min(weights(1:Np));
min_w_c=min(weights(Np+1:end));

            if print_type > 1
               full(inv(W)*A)
            end

return
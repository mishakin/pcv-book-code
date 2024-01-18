%% robustly estimate correction to dem, flatness
%
% [A,dx,v,W]=...
%   estimate_dem_slope_dependent...
%   (N,U,Np,Nr,Mc,poi,xa,weights,sigma_k,iter,out_print,type_robust);
% 
% N         = number of observations (points+second derivatives)
% U         = number of unknowns (dem)
% Np        = number of points
% Nr        = number of rows
% Mc        = number oc columns 
% poi       = Npx4 matrix of points (x,y,h,sigma), x,y, referring to grid
% xa        = approximate values for dem
% weights   = Nx1 vector of weights in [0,1]
% sigma_s   = standard deviation of first derivatives
% out_print = 0,1,2
% type_robust = 0,1,2,3, (00,01,10,11) for points and dem
%
% A         = NxU sparse designa matrix
% weights   = new weights
% dx        = delta dem
% v         = normalized residuals
% W         = weigths of obs. combining uncertainty and robust factor
%
% wf 10/2014
%

function [A,weights,dx,v,W]=...
    estimate_dem_robust_flat...
    (N,U,Np,Nr,Mc,poi,xa,weights,sigma_s,iter,type_robust,L1,...
    print_type,plot_type)
    

%% update weights indicating outliers
k=3;

if print_type > 0
    disp('estimate_dem_robust_flat')
    iter_typerobust_minW=[iter,type_robust,min(weights(Np+1,end))]
end
%% start estimation
% initiate size of A
if print_type > 0
    tic
end
A = spalloc(N,U,6*N);
b = zeros(N,1);
W = spalloc(N,N,N);

iter_switch=10;

% observations, normalized with 1/sqrt of variances
for p = 1:Np
    x = poi(p,1);
    y = poi(p,2);
    h = poi(p,3);
    w = sqrt(weights(p));
    s = poi(p,4);
    i = floor(x);
    j = floor(y);
    aa = x-i;
    bb = y-j;
    %delta_x
    A(p,i+j*Nr+1)       = (1-aa)*(1-bb)/s*w;
    A(p,i+j*Nr+2)       =    aa *(1-bb)/s*w;

    A(p,i+(j+1)*Nr+1)   = (1-aa)* bb   /s*w;
    A(p,i+(j+1)*Nr+2)   =    aa * bb   /s*w;
    delta_l =  h - ...
        ((1-aa)*(1-bb)*xa(i+j*Nr+1)    +...
        aa *(1-bb)*xa(i+j*Nr+2)    +...
        (1-aa)* bb   *xa(i+(j+1)*Nr+1)+...
        aa * bb   *xa(i+(j+1)*Nr+2));
    b(p)                =   delta_l   /s*w;
    W(p,p)=1/s*w;
end


if print_type > 0
    time_for_building_An=toc
end

%% priors
if print_type > 0
    tic
end
Icc = Np;                    % first index for column curvatures
Irr = Np + (Nr-1)*Mc;        % first index for row curvature
% coefficients for column curvature
% cc  1
%    -1
Jcc = [0,+1];  % index vector
Vcc = [-1 1];                    % coefficients
wcc = 1/sigma_s/sqrt(2);  % /norm(Vcc)
% coefficients for row curvature
% rr  1 -1
Jrr = [0,+Nr];  % index vector
Vrr = [-1  1];                    % coefficients
wrr = 1/sigma_s/sqrt(2); % /norm(Vrr)

% built up of priors in A
% column slope
for j=1:Mc
    for i=1:Nr-1
        IU  = i+(j-1)*Nr;
        %[i,j,IU]
        Icc = Icc+1;
        A(Icc,IU+Jcc)= Vcc * wcc * sqrt(weights(Icc));
        delta_cc     = 0-xa(IU+Jcc)' * Vcc';
        b(Icc)       = delta_cc * wcc * sqrt(weights(Icc));
        W(Icc,Icc)   = wcc * sqrt(weights(Icc));
    end
end
% row slope
for j=1:Mc-1
    for i=1:Nr
        IU  = i+(j-1)*Nr;
        %[i,j,IU]
        Irr = Irr+1;
        A(Irr,IU+Jrr)= Vrr * wrr * sqrt(weights(Irr));
        delta_rr     = 0-xa(IU+Jrr)' * Vrr';
        b(Irr)       = delta_rr * wrr * sqrt(weights(Irr));
        W(Irr,Irr)=wrr * sqrt(weights(Irr));
    end
end
%
if print_type > 0 
    time_for_building_A=toc
end
if print_type > 1
    if U < 1000
        Ab = 2*[full(A)]
        null_A_prior = U-rank(full(A(U+1:end,:)))
        rankA= rank(full(A(U+1:end,:)))
    end


        figure
        spy(A);


end
if print_type > 0
    time_for_building_Ae=toc
    start_time = cputime
end
%% solve
% no sort
% tic
% [C,R]=qr(A,b,0);
% dx = R\C;
% time_for_QR_original=toc
% % dx_unsorted = dx;


% sort
if print_type > 0
    tic
end
Nmatrix  = A'*A;

if print_type > 0
    size_N_matrix = size(Nmatrix)
    non_zeros_N_matrix =nnz(Nmatrix)
end
if iter==1 && plot_type > 0
    figure
    subplot(1,3,1)
    spy(Nmatrix)
end

%permx    = symrcm(Nmatrix);
permx    = symamd(Nmatrix);
Aperm    = A(:,permx);

if iter==1 && plot_type > 0
    subplot(1,3,2)
    spy(Aperm'*Aperm)
end

[C,R]    = qr(Aperm,b,0);
dxperm   = R\C;
Pm       = speye(U);
Pm       = Pm(:,permx);
dx       = Pm*dxperm;
if print_type > 0
    time_for_QR_sorted=toc
end
if iter==1 && plot_type > 0
    subplot(1,3,3)
    spy(R)
    fill_in = nnz(R)-(nnz(Nmatrix)/2+U/2);
    percentage_nnz   = nnz(R)/(U*(U+1)/2)
    relative_fill_in = fill_in/(nnz(Nmatrix)/2+U/2)-1
end

dx_sorted = dx;

% comparison
% comparison=[permx;dx_unsorted';dx_sorted';dx_sorted'-dx_unsorted';dxperm']

if U < 1600 && plot_type > 0
    Aqt = R\A';
    figure
    spy(Aqt')
end
% xo=xa;
% [ xa, err, iter, flag ] = jacobi_iterative_solution(A'*A, xa, A'*b, 200, 10^(-4))
% dx=xa-xo;

if print_type > 0
    time_for_solution=cputime-start_time
end


% residuals (normalized)
v = A*dx-b;
% estimated variance factor
if print_type > 0
    eso = norm(v)/sqrt(N-U)
end


eso_p=median(abs(v(1:Np)))*1.48;
eso_s=median(abs(v(Np+1:end)))*1.48;
if print_type > 0
    est_s0_s0p_s0k=[eso,eso_p,eso_s]
end

check_AtV_null=norm((A'*v));
max_v_points = max(abs(v(1:Np)));
max_v_slope  = max(abs(v(Np+1,end)));

% update weights indicating outliers

kp = k*eso_p;
kp = k;
ks = k*eso_s;
if print_type > 0
    criticalvalue_kp_sigmap_ks_sigmas=[kp,poi(1,4),ks,sigma_s]

    disp('vor Regewichtung')
    initial=[iter,type_robust,iter_switch,max(abs(v(Np+1:end)))]
end

if iter == 1 && type_robust > 0
    weights=ones(N,1);
else
    if type_robust > 0 && L1
        if type_robust == 1 | type_robust == 3
            weights(Np+1:end) = min(1,1./(abs(v(Np+1:end))/ks+eps));
            check=[iter,type_robust,min(weights(Np+1:end))];
        end
        if type_robust == 2 | type_robust ==3
            weights(1:Np) = min(1,1./(abs(v(1:Np))/kp+eps));
        end
    else
        if type_robust == 1 | type_robust == 3
            weights(Np+1:end) = exp(-abs(v(Np+1:end))/ks);
        end
        if type_robust == 2 | type_robust ==3
            weights(1:Np) = exp(-abs(v(1:Np))/kp);
        end
    end
end
abs_dx = norm(dx);
min_w_p=min(weights(1:Np));
min_w_c=min(weights(Np+1:end));

if print_type > 1
    full(inv(W)*A)
end

return
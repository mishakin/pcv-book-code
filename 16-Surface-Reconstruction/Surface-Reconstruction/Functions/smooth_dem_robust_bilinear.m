%% smooth_dem_robust_bilinear
%
% reconstruct surface with 
%   energy = sum_m rho(l_m-a_m)/sigma_n) + 
%             1/sigma_d sum_ij d_ii^2+2d_ij^2+djj^2
% l_m = observed heights, interpolated bilinearly
% a_m unknown grid heights, d_ij second derivatives
% 
% [dem_smoothed,S,Sigma,Np,Nr,Mc,v,A,weights,weights_f,W] = ...
%     smooth_dem_robust_bilinear...
%     (points,BB,delta_x,sigma_k,out_C,type_robust,...
%     out_in,print_type,plot_type)
%
% points      = Nx4-matrix of (x,y,h,sigma) values, 
% BB          = [xmin,ymin,xmax,ymax] bounding box of square grid
% delta_x     = grid size, squared 
%               xmax, ymax will be adapted, not decreased
%               points not inside BB will be discarded
% sigms_k     = std. dev. of curvature of surface, referring to delta_x
% out_C       = 0 no covariance matrix as output, (default if U>1600)
%               1 covariance matrix as output
% type_robust = 0,1,2 (00,01,10) for points and dem
% out_in      = boolean N-vector indicating inliers
% print_type  = 0 no output
%               1 timing
%               2 little output
%               2 much output
% plot_type   = 0,1,2,3
%
% dem_smoothed = Nr x Mc matrix of smoothed heigths
% S            = corresponding standard deviations, if required
% Sigma        = full covariance matrix, if required
% Np           = number of points
% Nr, Mc       = number of grid points in x- and y-direction
% v            = residuals (Np + (Nr-2)*Mc + Nr*(Mc-2) + (Nr-1)*(Mc-1)) for
%                points
%                second derivatives cc
%                second derivatives rr
%                torsion
% A            = design matrix (sparse) 
% weights      = weight-factors after robust estimation in [0,1]
% weights_f    = weight-factors after robust estimation (0 or 1)
% W            = weight of observations including both, 
%                 uncertainty and
% robust factor 
%
% grid is matrix Nr x Mc and stored column wise
%
% wf 7/2014, 4/2018
%
 
function [dem_smoothed,S,Sigma,Np,Nr,Mc,v,A,weights,weights_f,W] = ...
    smooth_dem_robust_bilinear...
    (points,BB,delta_x,sigma_k,out_C,type_robust,...
    out_in_0,print_type,plot_type)

    tic

%% number of points, preliminary --> final (which are in BB)
[Npp,tmp] = size(points);

% eliminate all points outside BB
xmin = BB(1);
ymin = BB(2);
xmax = BB(3);
ymax = BB(4);
Nr =ceil((xmax-xmin)/delta_x)+1;
Mc =ceil((ymax-ymin)/delta_x)+1;


ind = find((points(:,1)>=BB(1))&(points(:,1)<BB(3))...
          &(points(:,2)>=BB(2))&(points(:,2)<BB(4)));
Np  = size(ind,1);
poi = [(points(ind,1:2)-ones(Np,1)*[xmin,ymin])/delta_x,...
             points(ind,3:4)];
out_in=out_in_0(ind);


%% number of unknowns and observations incl. prior
U = Nr*Mc;
N = Np + (Nr-2)*Mc + Nr*(Mc-2) + (Nr-1)*(Mc-1);

%% Provide Meta data

if print_type>0 
    display(['Number of points              ',num2str(Npp)])
    if Npp ~= Np
        display(['Number of points reduced from ',num2str(Npp),' to ',num2str(Np)])
    end    
    display(['Grid size                     ',num2str(Nr),' x ',num2str(Mc)])   
    display(['Number of unknowns            ',num2str(U)])
    display(['Number of observations        ',num2str(N)])
    if type_robust(1) > 0        
        display(['Robust estimation        '])
    else
        display(['Non-robust estimation        '])
    end
    
end
%% do not give covariance matrix for U>40401.
            if U > 40401
                out_C = 0;
            end
            if print_type > 0 
                tic
            end
%% estimate iteratively

iter_switch=3;
if type_robust == 0
       N_iter = 1;
else
       N_iter = 2*iter_switch; % 3 non-robust, 3 robust
end

xa        = zeros(U,1);               % initial unknown parameters
weights   = ones(N,1);                % weight-factors for robust estimation
weights_f = weights;
sigma_k0  = 1*sigma_k;
gain=(sigma_k/sigma_k0)^(1/(N_iter-1));

time_for_checking=toc;

startcputime=cputime;

x_prev = zeros(U,1);

for iter = 1:N_iter
    if print_type > 0
        display(['Iteration: ',num2str(iter),' ------------------- '])
    end
    xa = zeros(U,1);
    L1 = 0;

if type_robust > 0 && iter <= iter_switch
    % the first three iterations use L1 (smoothed)
    L1=1;
else
    % the last three iterations use the exponential reweighting
    L1=0;
end
    
    if type_robust > 0
        sigma_k=sigma_k0*gain^(iter);     % adapt prior sigma
    end
    
    if iter == 1
        % store design matrix A and right hand sides b
        [A,b,permx,weights,xa,v,W,Nmatrix]=estimate_dem_robust_1...
            (N,U,Np,Nr,Mc,poi,weights,sigma_k,type_robust,...
            out_in,print_type,plot_type);
           
    else       
        % use design matrix A and right hand sides b
        [weights,xa,v,W,Nmatrix]=estimate_dem_robust_2_plus...
            (A,b,permx,N,U,Np,Nr,Mc,poi,weights,sigma_k,...
            type_robust,L1,out_in,print_type);
           
           
    end
    
     
    if print_type > 1
        iter_xa = [iter,norm(xa-x_prev)]
    end
    %% plot 
                if plot_type > 0
                    if iter == 1
                        if type_robust > 0
                            figure
                            hold on
                            X=([4*Nr:4*Nr+Nr-1]'*ones(1,Mc))*delta_x+BB(1);
                            Y=(ones(Nr,1)*[0:Mc-1])*delta_x+BB(2);
                            
                            imagesc(reshape(xa,Nr,Mc));
                            colormap(gray)
                            view(0,90)
                            axis equal;axis off
                            title('1.iteration, fitted dem as image','FontSize',16)
                          
                            xlim([0,70])
                            ylim([0,70])
                        end
                    else
                        figure
                        subplot(2,2,1)
                        hold on
                        X=([0:Nr-1]'*ones(1,Mc))*delta_x+BB(1);
                        Y=(ones(Nr,1)*[0:Mc-1])*delta_x+BB(2);
                        colormap(gray)
                        mesh(X,Y,reshape(xa,Nr,Mc),'EdgeColor',[0,0,0])
                        view([-29,65])
                        zlim([0,60])
                        grid off

                        title(strcat('iteration=',num2str(iter),' fitted dem'))

                        Icc = Np;                    % first index for column curvatures
                        Irr = Np + (Nr-2)*Mc;        % first index for row curvature
                        Irc = Np + (Nr-2)*Mc + Nr*(Mc-2);    % first index for torsion
                        subplot(2,2,2)
                        imagesc(reshape(weights(Icc+1:Irr),Nr-2,Mc));
                        % figure transposed ->
                        title(strcat('iteration=',num2str(iter),' weights z_r_r'))
                        
                        subplot(2,2,3)
                        if type_robust ~= 2
                            imagesc(reshape(weights(Irr+1:Irc),Nr,Mc-2));
                        else
                            rr=floor(sqrt(Np));
                            imagesc(reshape(1-weights(1:rr^2),rr,rr))
                        end
                        % figure transposed ->
                        title(strcat('iteration=',num2str(iter),' weights z_r_r'))
                        
                        subplot(2,2,4)
                        if type_robust ~= 2
                            imagesc(reshape(weights(Irc+1:end),Nr-1,Mc-1));
                        else
                            rr=floor(sqrt(Np));
                            imagesc(reshape(out_in(1:rr^2),rr,rr))
                        end
                        title(strcat('iteration=',num2str(iter),' weights z_r_c'))
                    end
                end
    x_prev = xa;
end
if type_robust > 0
    %final estimate
    if print_type > 0
        display('Last iteration -----------------')
    end
    xa        = zeros(U,1);
    weights_f = 0.9999*ceil(weights-0.5)+0.0001;
    if print_type > 2
        weights_final = weights_f(1:20)' 
    end
    xa=zeros(length(xa),1);
    
    % use design matrix A and right hand sides b
    [weights_f,dx,v,W,Nmatrix]=estimate_dem_robust_2_plus...
        (A,b,permx,N,U,Np,Nr,Mc,poi,weights_f,sigma_k,0,0,...
        out_in,print_type);
       
    %xa = xa + dx;
    xa = dx;
                if print_type > 1
                    iter_dx = [iter+1,norm(dx)] 
                end

                if plot_type > 0
                figure
                    subplot(2,2,1)
                    hold on
                    X=([0:Nr-1]'*ones(1,Mc))*delta_x+BB(1);
                    Y=(ones(Nr,1)*[0:Mc-1])*delta_x+BB(2);
                    colormap(gray)
                    mesh(X,Y,reshape(xa,Nr,Mc),'EdgeColor',[0,0,0])%  
                    %axis equal
                    view([-29,65])
                    zlim([0,60])
                    grid off
                    %alpha(0.3);
                    %mesh(ds);
                    %plot3(points(:,1),points(:,2),points(:,3),'.r','Marker
                    %Size',5);

                    title(strcat('last iteration, fitted dem'))

                    Icc = Np;                    % first index for column curvatures
                    Irr = Np + (Nr-2)*Mc;        % first index for row curvature
                    Irc = Np + (Nr-2)*Mc + Nr*(Mc-2);    % first index for torsion
                    subplot(2,2,2)
                    imagesc(reshape(weights_f(Icc+1:Irr),Nr-2,Mc));
                    % figure transposed ->
                    title(strcat('last iteration, weights z_c_c'))
                        
                    subplot(2,2,3)
                    if type_robust ~= 2
                        imagesc(reshape(weights_f(Irr+1:Irc),Nr,Mc-2));
                    else
                       rr=floor(sqrt(Np));
                       imagesc(reshape(1-weights_f(1:rr^2),rr,rr))
                    end
                    % figure transposed ->
                    title(strcat('last iteration, weights z_r_r'))
                        
                    subplot(2,2,4)       
                    if type_robust ~= 2
                       imagesc(reshape(weights_f(Irc+1:end),Nr-1,Mc-1));
                    else
                       rr=floor(sqrt(Np));
                       imagesc(reshape(out_in(1:rr^2),rr,rr))
                    end
                    title(strcat('last iteration, weights z_r_c'))
                end
end

dem_smoothed = reshape(xa,Nr,Mc);

%% determine covariance matrix for quality analysis
Sigma=0;
S = 0;

if out_C > 0    
                  
    Sigma = sparseinv(Nmatrix);     
                
    S     = full(reshape(diag(Sigma),Nr,Mc));
end
if print_type > 0
    display(['Total CPU time     = ',num2str(cputime-startcputime)])
end
return
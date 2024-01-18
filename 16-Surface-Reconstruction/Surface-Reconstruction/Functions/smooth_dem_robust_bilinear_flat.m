%% smooth_dem_robust_bilinear, sparse data on grid, flatness
%
% 
% function [dem_smoothed,S,Sigma,Np,Nr,Mc,v,A,weights,W] = ...
%    smooth_dem_robust_bilinear_flat...
%    (points,BB,delta_x,sigma_k,out_C,out_print,type_robust,type_data)
%
% points  = Nx4-matrix of (x,y,h,sigma) values, 
% BB      = [xmin,ymin,xmax,ymax] boundin,g box of square grid
% delta_x = grid size, squared 
%           xmax, ymax will be adapted, not decreased
%           points not incide BB will be discarded
% sigms_s = std. dev. of slope of surface, referring to delta_x
% out_C   = 0 no covariance matrix as output, (default if U>1600)
%           1 covariance matrix as output
% out_print = 0 no output
%             1 little output
%             2 much output
% type_robust = 0,1,2,3 (00,01,10,11) for points and dem
%
% dem_smoothe = Nr x Mc matrix of smoothed heigths
% S           = corresponding standard deviations, if required
% Sigma       = full covariance matrix. if required
% Np          = number of points
% Nr, Mc      = number of grid points in x- and y-direction
% v           = residuals (Np + (Nr-2)*Mc + Nr*(Mc-2)  for
%               points
%               first derivatives c
%               first derivatives r
% A           = design matrix (sparse) 
% weights     = weight-factors after robust estimation (0 or 1)
% W           = weight of observations including both, uncertainty and robust factor 
%
% grid is matrix Nr x Mc and stored column wise
%
% wf 10/2014
%

function [dem_smoothed,S,Sigma,Np,Nr,Mc,v,A,weights,W] = ...
    smooth_dem_robust_bilinear_flat...
    (points,BB,delta_x,sigma_s,out_C,type_robust,...
    print_type,plot_type)

% number of points, preliminary
[Npp,dim] = size(points);

if print_type > 0
    number_of_points_given = Npp
end


xmin = BB(1);
ymin = BB(2);
xmax = BB(3);
ymax = BB(4);
Nr =ceil((xmax-xmin)/delta_x)+1;
Mc =ceil((ymax-ymin)/delta_x)+1;


Np = 0;
for p=1:Npp
    if points(p,1) >= xmin && ...
            points(p,1) < xmax && ...
            points(p,2) >= ymin && ...
            points(p,2) < ymax
        Np=Np+1;
        poi(Np,:) = [(points(p,1)-xmin)/delta_x,(points(p,2)-ymin)/delta_x,points(p,3:4)];
    end
end

if print_type > 1
    number_of_points_in_BB = Np
end


% number of unknowns and observations incl. prior

U = Nr*Mc;
N = Np + (Nr-1)*Mc + Nr*(Mc-1);
if print_type > 0
    number_of_unknowns = U
end

% do not give covariance matrix for U>1600.
if U > 1600
    out_C = 0;
end
if print_type > 0
    tic
end
%% estimate iteratively

iter_switch=3;
if type_robust == 0
    N_iter = 1;
    first = 1;
else
    N_iter=6;
    first=1;
end
xa = zeros(U,1);
weights  = ones(N,1);
sigma_s0=sigma_s;
for iter = 1:N_iter
    if print_type > 0
        disp('Call estimate_dem_robust_flat')
        iter_typerobust=[iter,type_robust]
    end
    L1=0;
    if type_robust >0 && iter <= iter_switch
        L1=1;
        sigma_s=sigma_s0*4;
    else
        L1=0;
        sigma_s=sigma_s0;
    end
    sigma_s=sigma_s0*1.4^(N_iter-iter);
    [A,weights,dx,v,W]=estimate_dem_robust_flat...
        (N,U,Np,Nr,Mc,poi,xa,weights,sigma_s,iter,type_robust,L1,...
        print_type,plot_type);
    first=0;
    xa = xa + dx;
    

    if print_type > 0
        iter_dx = [iter,norm(dx)]
    end
    %%
    if  plot_type > 0
        if iter == 1
            figure(20)
            hold on
            X=([4*Nr:4*Nr+Nr-1]'*ones(1,Mc))*delta_x+BB(1);
            Y=(ones(Nr,1)*[0:Mc-1])*delta_x+BB(2);

            if type_data ~=5 && type_robust > 0
                colormap(gray)
                surf(X,Y,reshape(xa,Nr,Mc));
                alpha(0.3);
                %mesh(ds);

                plot3(points(:,1)+4*Nr*delta_x,points(:,2),points(:,3),'.r','MarkerSize',20);

            else
                subplot(1,3,2)
                %imagesc(reshape(xa,Nr,Mc));
                imshow(reshape(xa,Nr,Mc));
                colormap(gray)
                view(0,90)
                title('1.iteration')
            end
            title(strcat('iteration=',num2str(iter),' fitted dem'))

            axis equal
        else
            figure
            hold on
            X=([0:Nr-1]'*ones(1,Mc))*delta_x+BB(1);
            Y=(ones(Nr,1)*[0:Mc-1])*delta_x+BB(2);
            colormap(gray)
            surf(X,Y,reshape(xa,Nr,Mc));
            alpha(0.3);
            %mesh(ds);
            plot3(points(:,1),points(:,2),points(:,3),'.r','MarkerSize',20);

            title(strcat('iteration=',num2str(iter),' fitted dem'))

            if type_data ~= 5
                axis equal
            end
        end
    end

end
% if type_robust > 0
%     %final estimate
%     weights = 0.9999*ceil(weights-0.5)+0.0001;
%     weights_final = weights(1:20)'
%     xa=zeros(length(xa),1);
%     [A,weights,dx,v]=estimate_dem_robust_flat(N,U,Np,Nr,Mc,poi,xa,weights,sigma_s,1,out_print,0);
%     xa = xa + dx;
%     iter_dx = [iter+1,norm(dx)]
% end

dem_smoothed = reshape(xa,Nr,Mc);

%% determine covariance matrix
Sigma=0;
S = 0;
if print_type > 0
    tic
end
if out_C > 0
    Sigma = full(inv(A'*A));
    S     = full(reshape(diag(Sigma),Nr,Mc));
end
if print_type > 0
    time_for_inverse=toc
end


return
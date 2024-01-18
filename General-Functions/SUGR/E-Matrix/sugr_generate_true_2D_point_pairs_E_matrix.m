%% Generate true 2D point pairs from Relative Orientation,
% i.e., for given b, R generate N point pairs
%
% PP = sugr_generate_true_2D_point_pairs_E_matrix(b,R,N,boolean_r,Z0);
%
% b        = 3 x 1 vector 
% R        = 3x3 matrix R
% N         = number of points (in square/cube [-1,1]^2  / 3
% boolean_r = boolean: points should sit random
%                else: points sit in a square +- random (N should be square)
% Z0        = Z-shift of points from origin
%
% PP = point pairs
%      PP.h = N x 6 matrix of pairs of homogeneous point coordinates
%      PP.Crr = N x 4 x 4 sith 4 x 4 reduced CovM of point pairs
%      PP.type = 8 * ones(N,1)
% X  = true 3D coordinates
%
% Wolfgang Förstner 09/2011
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de
%
% See also sugr_E_Matrix

function [PP,X] = sugr_generate_true_2D_point_pairs_E_matrix(b,R,N,br,Z0)

K   = eye(3);
P1  = K *     [eye(3),zeros(3,1)];
P2  = K * R * [eye(3),-b];
% random points in unit square shifted by Z0 
X = zeros(N,4);
if br
    for n = 1:N
        %random points
        X(n,:)   = [rand(3,1)*2-1+[0,0,Z0]';1];
        % true image coordinates
        xs  = P1 * X(n,:)';
        xss = P2 * X(n,:)';
        xs  =xs/norm(xs);
        xss =xss/norm(xss);
        % struct of true point pairs      
        PP.h(n,:)     = [xs',xss'];
        PP.Crr(n,:,:) = zeros(4);
        PP.type(n)    = 8;
    end
    cov(X);
else
    M = ceil(sqrt(N));
    n = 0;
    for m = 1:M
        for k=1:M
           n = n+1;
           if n <= N
               X(n,:) = [-(M-1)/2+2*(m-1)/(M-1);...
                    -(M-1)/2+2*(k-1)/(M-1);...
                    Z0+rand(1)-0.5;...
                    1];
               % true image coordinates
                xs  = P1 * X(n,:)';
                xss = P2 * X(n,:)';
                xs  = xs/norm(xs);
                xss = xss/norm(xss);
                % struct of true point pairs      
                PP.h(n,:)     = [xs',xss'];
                PP.Crr(n,:,:) = zeros(4);
                PP.type(n)    = 8;
           end
        end
    end
end

figure('Color','w')
hold on
scatter3(X(:,1),X(:,2),X(:,3),'.r')
scatter3(0,0,0,'ob')
scatter3(b(1),b(2),b(3),'xb')
legend('object points','camera 1', 'camera 2')
xlabel('X')
ylabel('Y')
zlabel('Z')

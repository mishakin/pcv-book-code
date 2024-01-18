%% generate_points_on_plane(r,A);
%
% N     = number of points
% A     = plane (assumed to be not vertical, A3\=0))
% sigma = standard deviation (in all coordinates)
% r     = [xmin,xmax,ymin,ymax]
% RL    = 0 left plane
%       = 1 right plane
% rPCV  = 0 generate points
%       = 1 take coordiantes from book
%
% X      = list of noisy coordinates (structs)
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de

function X = generate_points_on_plane(N,A,sigma,r,RL,rPCV)

global plot_option

% Margin
m = 0.15;
% AX+BY+CZ+D=0;
Xlist = (1-2*m)*rand(N,1)*(r(2)-r(1))+r(1)+m*(r(2)-r(1));
Ylist = (1-2*m)*rand(N,1)*(r(4)-r(3))+r(3)+m*(r(4)-r(3));
Zlist = -(A(1)*Xlist + A(2)*Ylist + A(4)) / A(3);

X_true = [Xlist,Ylist,Zlist];

%% redundancy numbers for non-orthogonal fit
% 1-1/N-xq/sumxq-yq/sumyq

meanX = mean(Xlist);
meanY = mean(Ylist);
Xq = Xlist-meanX;
Yq = Ylist-meanY;
C = (N-1)*cov([Xq,Yq]);
D = [Xq,Yq]*inv(C)*[Xq,Yq]';                                               %#ok<MINV>
ri = ones(N,1)*(1-1/N)- diag(D);
deltai = 4./sqrt(ri);
sum(ri);
disp(['Redundancy numbers  : ',num2str(ri')])
disp(['Sensitivity factors : ',num2str(deltai')])

%% noise
if rPCV == 0
    Xe = X_true + rand_gauss([0,0,0]',sigma^2*eye(3),N)';
else    
    XLRn=[-0.8749 3.8996 2.3502 0.5370 2.1825 2.6343;...
        -1.6218 3.2287 2.0395 0.8666 1.3679 2.2913;...
        -1.6059 1.6108 1.8873 0.9658 3.6694 2.4097;...
        -1.4275 0.0609 2.0249 1.0134 2.3504 2.2827;...
        -1.2100 0.8542 2.4765 -0.0517 2.4755 2.9952;...
        0 0 0                 1.6559 -0.0646 2.0355];
    
    if RL==0
        Xe = XLRn(1:5,1:3);
    else
        Xe = XLRn(1:6,4:6);
    end
end

X.h   = zeros(N,4);
X.Crr = zeros(N,3,3);
for n=1:N
    Xc = sugr_Point_3D(Xe(n,:)',sigma^2*eye(3));
    X.h(n,:)     = Xc.h';
    X.Crr(n,:,:) = Xc.Crr;
end

if plot_option > 0
   z = -(A(1)*r(1:2) + A(2)*r(3:4) + A(4)) / A(3);
   
   for n=1:N
       plot3(Xe(n,1),Xe(n,2),Xe(n,3),'.k','MarkerSize',15);
   end
   plot3([r(1),r(2),r(2),r(1),r(1)],[r(3),r(3),r(4),r(4),r(3)],[z(1),z(2),z(2),z(1),z(1)],'--k');
   xlabel('X')
   ylabel('Y')
   zlabel('Z')
   axis equal
end

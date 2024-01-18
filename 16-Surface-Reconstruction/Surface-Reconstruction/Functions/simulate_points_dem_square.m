%% generate a set of points for dem-interpolation

function [points,BB,dx,sigma_k,sigma_h,oi,dem]=simulate_points_dem_square


%% square box, all observed
N=5;
% bounding box
BB=[0,0,N-1,N-1];
% grid size
factor=1;
dx = 1/factor;

sigma_h = 1.0;
sigma_k = 4/factor^(3/2);
k=0;
for i=1:N
    for j=1:N
        k=k+1;
        points(k,1:4)=[i,j,1,sigma_h];
    end
end

%%

xmin = BB(1);
ymin = BB(2);
xmax = BB(3);
ymax = BB(4);
Nr = ceil((xmax-xmin)/dx)+1;
Mc = ceil((ymax-ymin)/dx)+1;
%points(:,3)=(points(:,3)-mean(points(:,3)))*sqrt(Nr)/4;
% interpolation kernel
sigma = min(Nr,Mc)/5;

dem = zeros(Nr,Mc);
oi=ones(N^2,1);


return
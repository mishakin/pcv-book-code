%% generate a set of points for dem-interpolation
%
% dx   = grid size

function [points,BB,dx,sigma_k,sigma_s,dem,out_in]=simulate_points_dem_0(dx)


%% Example PCV (do not change)
% bounding box
BB=[-3,-3,14.0,10.0];

sigma_h = 0.001;
sigma_k = 0.001;
sigma_s = 0.1;

points = [...
     4.0, 4.0, 5,sigma_h;...
     2.0, 2.0, 1,sigma_h;...
     2.0, 6.0, 1,sigma_h;...
     6.0, 2.0, 1,sigma_h;...
     6.0, 6.0, 1,sigma_h;...
    10.0, 2.0, 3,sigma_h;...
    10.0, 6.0, 3,sigma_h ...
    ];
points(:,1:2)=points(:,1:2);
out_in=ones(size(points,1),1);
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
for i=1:Nr
    for j=1:Mc
        dem(i,j)=sum(points(:,3).*...
            exp(-1/2*((i-1-(points(:,1)-xmin)/dx).^2+(j-1-(points(:,2)-ymin)/dx).^2)/sigma^2));
    end
end

return
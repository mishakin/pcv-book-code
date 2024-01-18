%% generate a set of points for dem-interpolation

function [points,BB,dx,sigma_k,sigma_h,out_in,dem]=simulate_points_dem_0


%% Example PCV (do not change)
% bounding box
BB=[0,0,9,6.000]*2/3;
% grid size
factor=1;
dx = 1/factor;

sigma_h = 1.0;
sigma_k = 4/factor^(3/2);

points = [...
    1.2, 1.5, 2,sigma_h;...
    1.8, 1.2, 3,sigma_h;...
    5.4, 2.4, 4,sigma_h;...
    2.1, 4.8, 5,sigma_h;...
    8.7, 3.3, 6,sigma_h;...
    6.0, 5.4, 1,sigma_h ...
    ];

points(:,1:2)=points(:,1:2)*2/3;

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


out_in=ones(6,1);

return
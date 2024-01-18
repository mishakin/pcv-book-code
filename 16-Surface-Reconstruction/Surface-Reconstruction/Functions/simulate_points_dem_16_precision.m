%% generate a set of points for dem-precision analysis
%
% N           = number of points
% sigma_h     = precision of points
% Resolution  = Number of grid cells in each coordinate axis
%
% points      = Nx4 matrix of points
%               [X,Y,Z,sigma]
% BB          = bounding box
% dx          = sampling distance
% dem         = intended dem (not used here)
% bin_im      = binary image of points
% dist_t      = distance transformation of points
%
% wf 8/2015


function [points,BB,dx,dem,bin_im,dist_t]=simulate_points_dem_16_precision(N,sigma_h,R)


%% Example PCV (do not change)
% bounding box
BB=[0,0,R,R];
% grid size
factor=1;
dx = 1/factor;

p    = (1+rand(N,2)*(R-2));
M=randn(4);
a =[ones(N,1) (2*p(:,1)-R)/R (4*p(:,1).^2-4*p(:,1)*R-R^2)/(2*R^2) ...
    (2*p(:,1)-R).*(2*p(:,1).^2-2*p(:,1)*R-R^2)/R^3];
b =[ones(N,1) (2*p(:,2)-R)/R (4*p(:,2).^2-4*p(:,2)*R-R^2)/(2*R^2) ...
    (2*p(:,2)-R).*(2*p(:,2).^2-2*p(:,2)*R-R^2)/R^3];
z = 2*diag(a*M*b');
points = [p,z,ones(N,1)*sigma_h];
 
k=0;
for i=[1,round(R/2),R-1]
    for j=[1,round(R/2),R-1]
        k=k+1;
        points(k,:)=[i,j,2*[1 (2*i-R)/(2*R) 1/2*(4*i^2-4*i*R-R^2)/R^2 ...
            (2*i-R)*(2*i^2-2*i*R-R^2)/R^3]*M*...
            [1 (2*j-R)/(2*R) 1/2*(4*j^2-4*j*R-R^2)/R^2 ...
            (2*j-R)*(2*j^2-2*j*R-R^2)/R^3]',sigma_h];
    end
end



%% plot points and distance transform

xmin = BB(1);
ymin = BB(2);
xmax = BB(3);
ymax = BB(4);

Nr = ceil((xmax-xmin)/dx)+1;
Mc = ceil((ymax-ymin)/dx)+1;

dem = zeros(Nr,Mc);

bin_im = dem;
for n=1:N
    i=ceil(points(n,1)*factor);
    j=ceil(points(n,2)*factor);
    bin_im(R+1-j,i)=1;
end

ss = plot_init;
figure('Color','w','Position',[0.1*ss(1),0.2*ss(2),0.7*ss(1),0.6*ss(2)])
subplot(1,2,1)
imshow(bin_im)
axis equal
dist_t=bwdist(bin_im,'euclidean');
subplot(1,2,2)
imagesc(sqrt(dist_t));
title('points and distance transform')
axis equal
axis off



return
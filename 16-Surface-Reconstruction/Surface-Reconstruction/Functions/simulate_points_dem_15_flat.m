%% generate a set of points for dem-interpolation

function [points,BB,dx,sigma_k,sigma_s,out_in,dem]=simulate_points_dem_15(d,n)

d=2*d
N=(n-1)*d+1;
%% 
% bounding box
BB=[0,0,N+1,N+1];
% grid size
dx = 1;

sigma_h = 0.03725;
sigma_k = 1;
sigma_s = 1;

points = zeros(n*n,4);
k=0;
for i=1:n
    for j=1:n
        k=k+1;
        points(k,:)=[1+(i-1)*d,1+(j-1)*d,0,sigma_h];
    end
end
points(:,1:2)=points(:,1:2);
Np=size(points,1);
out_in=ones(Np,1);

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

figure
hold on
X=([0:Nr-1]'*ones(1,Mc))*dx+BB(1);
Y=(ones(Nr,1)*[0:Mc-1])*dx+BB(2);
surf(X,Y,dem)
colormap(gray)
alpha(0.3)
%mesh(dem)

%plot3(points(:,1)/dx+1,points(:,2)/dx+1,points(:,3),'.r','MarkerSize',20)
plot3(points(:,1),points(:,2),points(:,3),'.r','MarkerSize',15)
axis equal
xlim([xmin-1,xmax+1]);
ylim([ymin-1,ymax+1]);
title(strcat('original grid, n=',num2str(n)));

return
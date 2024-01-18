%% simulate_points_dem_6: generate test image of Uwe 1994
%
% [points,BB,delta_x,sigma_k,sigma_s,sigma,out_in,dem] = ...
%        simulate_points_dem_6(N,g_min,delta_g,sigma,outlier_percentage)
%
% Input:
%    N:                   int, mesh size in x- and y-direction
%    g_min:               double, min high level
%    delta_g:             double, hight steps
%    sigma:               double, std noise
%    outlier_percentage:  double \in [0,1], fraction of outlies
%
% Output
%    points:    double (N*N)x4, coordinates of dem grid points, [i, j, height(i,j), std]
%    BB:        int 1x4, bounding box of dem
%    delta_x:   double, gridsize/spacing
%    sigma_k:   double, std curvature   
%    sigma:     double, std noise
%    out_in:    boolean (N*N)x1, outlier mask
%    dem:       double NxN, DEM
%
% Wolfgang Förstner 07/14
% last changes: Susanne Wenzel 09/16, wf 4/2018
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de
%
function [points,BB,delta_x,sigma_k,out_in,dem] = ...
        simulate_points_dem_6(N,g_min,delta_g,sigma,outlier_percentage)

%% plot settings
ss = plot_init;

BB=[1,1,N,N];

%% set sigma of curvatures
sigma_k = 1.0*sigma;

%% start simulation
A=ones(N,N)*g_min;

%% squares: center (16,16)
cent_x=round(2/7*N);
cent_y=round(2/7*N);
% lowest: half width = 12, z = 40
h = round(12/70*N);
z = g_min+delta_g;
for i=cent_x-h:cent_x+h
    for j=cent_y-h:cent_y+h
        A(i,j)=z;
    end
end
% lowest: half width = 8, z = 60
h = round(8/70*N);
z = g_min+2*g_min;
for i=cent_x-h:cent_x+h
    for j=cent_y-h:cent_y+h
        A(i,j)=z;
    end
end
% lowest: half width = 4, z = 80
h = round(4/70*N);
z = g_min+3*delta_g;
for i=cent_x-h:cent_x+h
    for j=cent_y-h:cent_y+h
        A(i,j)=z;
    end
end

%% cylinder, sphere, block, peak
cent_x = round(2/7*N);
cent_y = round(5/7*N);
r      = round(13/70*N);

%% sphere
for i = cent_x-r:cent_x
    for j = cent_y-r:cent_y
      dx = (i-cent_x)/r;
      dy = (j-cent_y)/r;
      dr = sqrt(dx^2+dy^2);
      if dr < 1;
          A(i,j) = g_min + 3*delta_g*(1-sqrt(dx^2+dy^2));
      end
    end
end

%% box with peak
for i = cent_x-r:cent_x
    for j = cent_y:cent_y+r
      if dr < 1;
          A(i,j)=g_min + 3*delta_g;
      end
    end
end
A(cent_x-round(r/2),cent_y+round(r/2)) = g_min+4*delta_g;

%% cylinder
for i = cent_x:cent_x+r
    for j = cent_y-r:cent_y+r
      dx = (i-cent_x)/r;
      if dr < 1;
          A(i,j) = g_min + 3*delta_g*sqrt(1-dx^2);
      end
    end
end

%% diagonal cross
cent_x = round(5/7*N);
cent_y = round(2/7*N);
r      = round(13/70*N);
for i = cent_x-r:cent_x+r
    for j = cent_y-r:cent_y+r
        di = (i-cent_x)/r;
        dj = (j-cent_y)/r;
      A(i,j) = max(A(i,j),g_min+min([(1-di-dj),(1+di+dj)])*3*delta_g);
      A(i,j) = max(A(i,j),g_min+min([(1-di+dj),(1+di-dj)])*3*delta_g);
    end
end

%% torus and box
cent_x = round(5/7*N);
cent_y = round(5/7*N);
b      = round(12/70*N);
R      = round(12/70*N);
r      = round(03/70*N);
rm     = (R+r)/2;
for i = cent_x-b:cent_x+b
    for j = cent_y-b:cent_y+b
      dx = (i-cent_x);
      dy = (j-cent_y);
      dr = sqrt(dx^2+dy^2);
      ddr = 2*(dr - rm)/(R-r);
      if dr >= r && dr <= R
          A(i,j) = g_min + 3*delta_g*(1-sqrt(ddr^2));
      end
    end
end
for i = cent_x-b:cent_x-rm
    for j = cent_x-b:cent_y+b
        A(i,j) = g_min+3*delta_g;
    end
end

%% add noise
M = (rand(N,N) < outlier_percentage);
A = A + randn(N,N)*sigma + M .*(rand(N,N)-0.5)*8*delta_g;

%% Check
figure('name','original dem as image','color','w',...
    'Position',[0.02*ss(1),0.02*ss(2),0.3*ss(1),0.4*ss(2)]);
imagesc(A)
axis equal;axis off
colormap(gray)
title('original dem as image','FontSize',16)

%% output
dem = A;
delta_x = 1;
k = 0;
points = zeros(BB(3)*BB(4),4);
out_in = ones(BB(3)*BB(4),1);
for i=BB(1):BB(3)
    for j=BB(2):BB(4)
        k=k+1;
        points(k,:) = [i,j,A(i,j),sigma];
        out_in(k) = 1-M(i,j);
    end
end
dem = dem(BB(1):BB(3),BB(2):BB(4));
dem = dem(BB(1):BB(3),BB(2):BB(4));

figure('name','original dem','color','w',...
    'Position',[0.02*ss(1),0.52*ss(2),0.3*ss(1),0.4*ss(2)]);
plot_surface(dem,BB,delta_x,'plotfun',@mesh,'view',[-65,29],'alpha',1);
axis equal
title('original dem as mesh','FontSize',16)
return
%% read ply-file
%
% [points,BB,delta_x,sigma_k]=simulate_points_dem_10(file_name);
%
% file_name     = name of ply-file
%                 given in '..\Data'
%                 fa2_aurelo_result_pyra0_ausgeschnitten-1.ply
%                        (used in Fig. 16.23)
%                 fa2_aurelo_result_pyra0_ausgeschnitten-4.ply
%                        (larger data set)
% Nbin          = resolution
%
% points        = N x 4 vector of points and sigma
% BB            = bounding box
% delta_x       = grid distance
% sigma_k       = sigma of curvatures
%
% Wolfgang Förstner 07/14
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

function [points,BB,delta_x,sigma_k,XYZ_ori]=simulate_points_dem_10(file_name,Nbin)

%% plot settings
ss = plot_init;

%% read file
pts = ply_read(file_name);

%% extract coordinates
XYZ_ori = [pts.vertex.x,pts.vertex.y,pts.vertex.z];
N = size(XYZ_ori,1);

%% possibly rotate points to see them from the side
% % determine principal axes
% C = cov(XYZ_ori);
% [R,D] = eigs(C);
% R = R*diag(diag(sign(R)));
% % normalize
% XYZ = XYZ_ori*R;
XYZ = XYZ_ori;

%%
xmin = min(XYZ(:,1));
xmax = max(XYZ(:,1));
ymin = min(XYZ(:,2));
ymax = max(XYZ(:,2));
frame_size = 1.1;
bord = (frame_size-1)/2;
Range_x = (xmax-xmin);
Range_y = (ymax-ymin);
delta_x = min(Range_x,Range_y)/(Nbin-1);
Nr = ceil(Range_x*frame_size/delta_x);
Mc = ceil(Range_y*frame_size/delta_x);
BB = [xmin-bord*Range_x,               ymin-bord*Range_y,...
      xmin-bord*Range_x+(Nr-1)*delta_x,ymin-bord*Range_y+(Mc-1)*delta_x];

% sigma_n
factor_n = 0.1;
sigma_n = factor_n*max([Range_x,Range_y])/sqrt(N);
sigma_k = 3.0 * sigma_n;

% prepare points
points = [XYZ(:,1),XYZ(:,2),XYZ(:,3),sigma_n*ones(N,1)];

%%
figure('name','Given points','color','w',...
    'Position',[0.1*ss(1),0.3*ss(2),0.3*ss(1),0.4*ss(2)]);
plot3(points(:,1),points(:,2),points(:,3),'.k')
axis equal;axis off;hold on
xlabel('x');ylabel('y');zlabel('z');
view([0,90])
title('Given points')
plot3([BB(1),BB(1),BB(3),BB(3),BB(1)],[BB(2),BB(4),BB(4),BB(2),BB(2)],[0,0,0,0,0],'-k');




return

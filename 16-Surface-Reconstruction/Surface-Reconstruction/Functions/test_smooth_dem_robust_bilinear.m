%% test dem smoothing

% clear all
clearvars
close all
% type_data: type of generated data
% 0     a set of given points
% 1     a roof parallel to the grid
% 2     a roof skew to the grid
% 3     a roof skew to the grid + a flat roof, points can be thinned out
% 4     same as 3, but with points not sitting on a grid
% 5     read an image
% 6     Uwe's test image
% 7     cake
% 8     1/4 cake
% 9     constant height with outliers
% 10    read point cloud
% 11    flat
% 12    cone
% 13    cake
% 14    flat with outliers
% 15    flat 4x4 points
%
% wf 7/2014

init_rand = 6;
type_data = 0;
type_robust = 0
%             0 not robust
%             1 only dem
%             2 only points
%             3 both

out_C = 0;
out_print = 0;
print_type = 0;
plot_type  = 0;
Tol =0.15;

init_rand_seed(init_rand);


%% generate dem point cloud
switch type_data
    case 0
        [points,BB,delta_x,sigma_k,sigma_h,dem]=simulate_points_dem_0
        out_in =ones(size(points),1);
    case 1
        [points,BB,delta_x,sigma_k,dem]=simulate_points_dem_1;
    case 2
        [points,BB,delta_x,sigma_k,dem]=simulate_points_dem_2;       
    case 3
        [points,BB,delta_x,sigma_k,dem]=simulate_points_dem_3;
    case 4
        [points,BB,delta_x,sigma_k,dem]=simulate_points_dem_4;
    case 5
        im_name = 'IMG_2092_4.JPG'; 
        im_name = 'IMG_2726-2.JPG';  
        im_name = 'IMG_8004-4.JPG';   
        % im_name = 'IMG_2425-2.JPG';  % Figur Paris  
        %im_name = 'IMG_5977-4.JPG'; % Zebra
        %im_name = 'IMG_7452-4.JPG'; % Katze
        im_name = 'IMG_8942-4.JPG'; % Fassade Prag
        [points,BB,delta_x,sigma_k,sigma_h,dem]=simulate_points_dem_5(im_name);
     case 6
        [points,BB,delta_x,sigma_k,sigma_s,sigma_h,out_in,dem]=simulate_points_dem_6(70,10,10,1,0.05);
        [points,BB,delta_x,sigma_k,sigma_s,sigma_h,out_in,dem]=simulate_points_dem_6(70,10,10,0.5,0.05);
    case 7
        [points,BB,delta_x,sigma_k,dem]=simulate_points_dem_7(30,15,15,1);
        dem0=dem;
     case 8
        [points,BB,delta_x,sigma_k,dem]=simulate_points_dem_8(10,15,15,1);
        dem0=dem;
     case 9
         sigma=0.001;
        [points,BB,delta_x,sigma_k,dem,out_in,dz]=simulate_points_dem_9(2000,0.15,sigma);
        dem0=dem;
     case 10
        [points,BB,delta_x,sigma_k,XYZo]=...
            simulate_points_dem_10('fa2_aurelo_result_pyra0_ausgeschnitten-3.ply',40);
        
    case 11
        [points,BB,delta_x,sigma_k,sigma_s,dem]=simulate_points_dem_0_flat
    case 12
        [points,BB,delta_x,sigma_k,sigma_s,dem]=simulate_points_dem_cone
    case 13
        [points,BB,delta_x,sigma_k,sigma_s,dem]=simulate_points_dem_steps(60,4,5)
    case 14
         sigma=0.001;
        [points,BB,delta_x,sigma_k,dem,out_in,dem]=simulate_points_dem_14(2000,0.15,sigma);
        dem0=dem;
    
    case 15
        d=8;
        n=4;
        [points,BB,delta_x,sigma_k,sigma_s,out_in,dem]=simulate_points_dem_15_flat(d,n)        
        out_C = 1;
end

Np = size(points,1);


%% interpolate

starttime = cputime
[ds,S,Sigma,Np,Nr,Mc,ver,A,w,w_f,W] =...
    smooth_dem_robust_bilinear...
    (points,BB,delta_x,sigma_k,out_C,type_robust,type_data,...
    out_in,print_type,plot_type)
complete_time_for_solution=cputime-starttime
    
% figure
%        sc = max(ds(:))-min(ds(:));
%        imshow((ds-min(ds(:)))/sc*0.95);
%        
% if out_C > 0
%     figure
%     f=max(S(:));
%     mesh(sqrt(S/f*min(Nr,Mc)/3)*10)
%     hold on
%     title('standard deviations')
%     if type_data ~= 5 
%         axis equal
%     end
% end
%%
figure
hold on
X=([0:Nr-1]'*ones(1,Mc))*delta_x+BB(1);
Y=(ones(Nr,1)*[0:Mc-1])*delta_x+BB(2);

if type_data ~= 5
   colormap(gray)
    shade = [-1 0; 0 1]/sqrt(2);
    shade = -[ 0 -1; 1 0]/sqrt(2);
    col=conv2(ds,shade,'same');
    %col=(1+col./sqrt(1+col.^2))/2;
    colo=min(1, max(0,tanh(100*(col-0)/(max(col(:))-min(col(:))))));
    surf(X,Y,ds,colo,'FaceLighting','gouraud','EdgeColor','none')
    alpha(0.3);
    shading interp 
    %mesh(X,Y,ds);
    
    plot3(points(:,1),points(:,2),points(:,3),'.r','MarkerSize',15)
    view=[-20,30]
else
    figure %subplot(1,3,3)
    %imagesc([dem,ds])
    imshow([dem,ds])
    colormap(gray)
end

title('fitted dem with given points')
       
if type_data ~= 5
    axis equal
end
%%
figure
if type_data == 0
    subplot(1,2,1)

    hold on
    surf(X,Y,ds); %,'EdgeColor',[0,0,0])% ,129*ones(size(S))
    colormap(gray)

    title('fitted dem - min gradient')
    alpha(0.3)
    col=zeros(Np,2);
    for n=1:Np
        plot3(points(n,1),points(n,2),points(n,3),'.k','MarkerSize',15)
    end
    view([-29,65])
    xlim([BB(1),BB(3)])
    ylim([BB(2),BB(4)])
    zlim([0,10])
    axis equal

    subplot(1,2,2)
    hold on
    surf(X,Y,sqrt(S)*sigma_h,'EdgeColor',[0,0,0])
    colormap(cool);
    title('fitted dem - sigmas')
    col=zeros(Np,2);
    alpha(0.3)
    for n=1:Np
        z=interpolate_bilinear(S,points(n,1),points(n,2),delta_x,BB,sigma_h);
        plot3(points(n,1),points(n,2),z,'.k','MarkerSize',15)
    end
    plot3(0,0,0,'.k')
    view([-29,65])
    xlim([BB(1),BB(3)])
    ylim([BB(2),BB(4)])
    zlim([-1,7])
    axis equal

end

%%

figure

mesh(ds,'EdgeColor',[0,0,0])
axis equal
        VIEW([-33,63])
        zlim([-30,70])
grid off
%%
% figure
% hold on
% mesh_edge_color(ds)
% view([-68,57])
% axis equal

%plot_edges_dem(Np, Nr, Mc, points, out_in, w, w_f, Tol);
% for p=1:Np
%     if out_in(p)
%         [p,points(p,1),points(p,2),w(p),w_f(p),abs(dz(p))/(3*sigma)];
%     end
% end


figure
hold on
mesh(S);
title('Heigh \sigma_z')
view([-33,63])
sigma_centre=S(1+(n-1)*d,1+(n-1)*d)
return 
%%
if type_data ~= 5 
        axis equal
end
figure
imagesc(ds-dem)
colormap(gray)

figure
mesh(ds-dem)
colormap(gray)

full(inv(W)*A);

figure
mesh(ds)
colormap(gray);
axis equal

number_of_outliers = sum(out_in)


return
%%
if type_data == 1
    Icc = Np;                    % first index for column curvatures
    Irr = Np + (Nr-2)*Mc;        % first index for row curvature
    Irc = Np + (Nr-2)*Mc + Nr*(Mc-2);    % first index for torsion
    % plot residuals
    figure
    % heights
    subplot(1,4,1)
    res=zeros(Nr,Mc);
    vg = reshape(ver(1:Np),Nr-2,Mc-2);
    for i=1:Nr-2
        for j=1:Mc-2
            res(i+1,j+1)=vg(i,j);
        end
    end
    mesh(-res(:)');
    title('residuals - heights')
    
    % curvatures cc
    subplot(1,4,2)
    res=zeros(Nr,Mc);
    vg = reshape(ver(Icc+1:Irr),Nr-2,Mc);
    for i=1:Nr-2
        for j=1:Mc
            res(i+1,j)=vg(i,j);
        end
    end
    mesh(-res'*sqrt(sigma_k));
    title('residuals - curvatures cc')
    
    % curvatures rr
    
    subplot(1,4,3)
    res=zeros(Nr,Mc);
    vg = reshape(ver(Irr+1:Irc),Nr,Mc-2);
    for i=1:Nr
        for j=1:Mc-2
            res(i,j+1)=vg(i,j);
        end
    end
    mesh(-res'*sqrt(sigma_k));
    title('residuals - curvatures rr')
    
    % torsion rc
    
    subplot(1,4,4)
    res=zeros(Nr,Mc);
    vg = reshape(ver(Irc+1:end),Nr-1,Mc-1);
    for i=1:Nr-1
        for j=1:Mc-1
            res(i,j)=vg(i,j);
        end
    end
    mesh(-res'*sqrt(sigma_k));
    title('residuals - torsion rc')
end

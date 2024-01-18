% demo_h_determine_homography_from_two_images: read/interactively determine homography
%
% [H,X,Cxx] = demo_h_determine_homography_from_two_images...
%    (tau,sigma,sidth,sigma_n,Image_l,Image_r,coordinatesigma,coordinate_name,...
%    ,magnification,readX);
%
% H homography from left to right
% X    = 4 x 4-matrix, rows control point coordinates
% Cxx  = 4 x 4 x 4 tensor with cov. matrices
%
% Wolfgang Förstner 06/18
% wfoerstn@uni-bonn.de

function [H,X] = determine_Homography_2D_from_two_images(...
    Image_l,Image_r,readX)

global coordinate_sigma_x
global coordinate_name

addpath(genpath('../../General-Functions'))
ss = plot_init;

f1 = figure('Color','w','Name', 'left image','Position',[50,100,0.45*ss(1),0.45*ss(1)]);
imshow(Image_l);
hold on
f2 = figure('Color','w','Name', 'right image','Position',[100+0.45*ss(1),100,0.45*ss(1),0.45*ss(1)]);
imshow(Image_r);
hold on

if readX == 0
    X = zeros(4,4);
    Cxx = zeros(4,4,4);
    for n=1:4
        figure(f1)
        hold on
        disp('measure point in left image (1)') 
        %imagesc(Image_l);
        %% wait until a button is pressed
        keydown = waitforbuttonpress;
        % if a keyboard button was pressed, stop cycle
        if keydown == 1
            break
        end
        %% grab mouse position
        x = get(gca,'currentPoint');
        % determina local accuracy
        z = x(1,1:2)';

        Covl = coordinate_sigma_x*eye(2);

        %% plot point in left image
        plot_square_with_background(z(1),z(2),10);
        hold off
        
        X(n,1:2)= x(1,1:2);

        figure(f2)
        hold on
        
        disp('measure corresponding point in right image (2)') 
        %% wait until a button is pressed
        keydown = waitforbuttonpress;
        % if a keyboard button was pressed, stop cycle
        if keydown == 1
            break
        end
        %% grab mouse position
        x = get(gca,'currentPoint');
        z = x(1,1:2)';
        
        Covr = coordinate_sigma_x*eye(2);

        %plot_square_with_background(x(1,1),x(1,2),10);
        plot_square_with_background(z(1),z(2),10);
        X(n,3:4) = x(1,1:2);

        Cxx(n,:,:)= [Covl zeros(2);zeros(2) Covr];
    end

    save(coordinate_name,'X','Cxx');

else
    load(coordinate_name,'X','Cxx');
end
% conditioning
Tl = sugr_Homography_2D(X(:,1:2));
Tr = sugr_Homography_2D(X(:,3:4));
for n = 1 : 4
    pl = sugr_Point_2D(X(n,1:2)', squeeze(Cxx(n,1:2,1:2))); 
    pl = sugr_transform_with_Homography_2D(Tl,pl);
    
    pr = sugr_Point_2D(X(n,3:4)', squeeze(Cxx(n,3:4,3:4))); 
    pr = sugr_transform_with_Homography_2D(Tr,pr);
    %
    PPc.h(n,:)     = [pl.h', pr.h'];
    PPc.Crr(n,:,:) = [pl.Crr zeros(2);zeros(2) pr.Crr]; 
end

% estimation of H (conditioned)
Ha = sugr_estimation_algebraic_Homography_2D_from_point_pairs(PPc);
T = 0.1;
maxiter = 3;
[Hc,~,~] = sugr_estimation_ml_Homography_2D_from_point_pairs(PPc,Ha,T,maxiter);

% uncondition xr = H xl, 
% Tr xr = (Tr H Tl^{-1} Tl) xl 
% H = inv(Tr) * H * Tl
H = sugr_transform_Homography_2D(Hc, inv(Tr.H), Tl.H);


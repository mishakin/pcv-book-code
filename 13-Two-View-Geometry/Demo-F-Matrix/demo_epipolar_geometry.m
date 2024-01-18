%% Epipolar geometry
%
% The script allows to determine the relative orientation (F-matrix) 
% of two images from 8 corresponding points and to show the 
% epipolar lines of points in the right image 
% togehter with their uncertainty.
%
% The uncertainty of the measured points can be chosen to be 
%  1. the rounding error(appr. 1/3 pixel) or 
%  2. derived from the image content. 
%     Then the script can be used to explore 
%     the uncertainty of matching based on the sstructure tensor of
%     the image, see PCV (13.93)
% This is controlled by the variable 'cov_im'

% The uncertainty of the epipoles (if visible) and 
% the epipolar lines is visualized. 
%
% The 8 corresponding points and the fundamental matrix may be stored, 
% and read in later. As soon as 8 correspoinding points are measured, 
% the old file will be owerwritten.
%
% 1. Choose two corresponding images 
%    (They are assumed to be in the folder ../Data/images)
%    Six image pairs are prepared, so only the number is to be chosen
%    (see the parameter 'image_pair')
% 2. Decide on whether you want to measure 8 corresponidng points or 
%    read data from file (see the parameter 'readX')
%    (they are lying in the folder /Data/ObservedImageCoordinates
% 3. If chosen: measure the corresponding points 
%    one after the other in the first and the second image
% 4. Measure individual points in the first image (now number 3) 
%    to get the epipolar line in the second image (number 4) visualized.
%    In addition the expected precision of the measured 
%    parallaxes/disparities px' and py'is shown in image 3, 
%    and the standard deviations of px' and py' are give on the console.
% 5. Step 4 can be repeated until right mouse button is pressed.
%
% The folder '../Data/images' contains the original images.
% The folder '../Data/ObservedImageCoordinates' contains the measured points
% for each image pair twice
%    - as 'XXX.mat'-file which can be overwritten by measuring the
%      8 corresponding points
%    - as 'XXX-original.mat'-file, which contains already measured points,
%      which can be used to recover from wrongly measured points.
%
% Wolfgang Förstner 4/2018
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 05/18
% wenzel@igg.uni-bonn.de

function demo_epipolar_geometry()

addpath(genpath('../../General-Functions'))
addpath('../Functions')
addpath(genpath('../Data'))

ss = plot_init();
close all

% magnification = magnification factor for plotting 
plot_params.magnification_m = 300; % factor for measured point
plot_params.magnification   =  3; % factor for epipolar geometry

%% Control parameters
% choose number of image pair
image_pair = 1;

% choose whether to measure 8 image pairs or read data from file
% if readX = 0 then measure image correspondencies
%          = 1 read them from file
readX = 1;

% choose what uncertainty is shown for the expected parallaxes/disparities
% cov_im = 0 the uncertainty is the rounding error
%        ~=0 the uncertainty is derived from the local structure tensor
cov_im = 1;

%%
switch image_pair
    case 1
        %% Cathedral baseline in viewing direction
        Image_name_l = 'IMG_9473.JPG';
        Image_name_r = 'IMG_9474.JPG';
        plot_params.magnification_m = 200; % factor for measured point
        plot_params.magnification = 4;   % factor for epipolar geometry
    case 2
        %% Church
        Image_name_l = '008.JPG';
        Image_name_r = '009.JPG';
        plot_params.magnification_m = 100; % factor for measured point
        plot_params.magnification = 1;  % factor for epipolar geometry
    case 3
        %% Upside down building
        Image_name_l = 'IMG_1936.JPG';
        Image_name_r = 'IMG_1937.JPG';
        plot_params.magnification_m = 300; % factor for measured point
        plot_params.magnification = 6;   % factor for epipolar geometry
    case 4
        %% Cathedral b/w
        Image_name_l = 'IMG_9473-bw.JPG';
        Image_name_r = 'IMG_9474-bw.JPG';
        plot_params.magnification_m = 200; % factor for measured point
        plot_params.magnification = 4;   % factor for epipolar geometry
    case 5
        %% corn rigth left
        Image_name_l = '1336-34.jpg';
        Image_name_r = '1336-95.jpg';
        plot_params.magnification_m = 200; % factor for measured point
        plot_params.magnification = 2;   % factor for epipolar geometry
    case 6
        %% corn up/down
        Image_name_l=  '1337-1.jpg';
        Image_name_r = '1337-2.jpg';
        plot_params.magnification_m = 200; % factor for measured point
        plot_params.magnification = 2;   % factor for epipolar geometry
    otherwise
        disp('wrong number of image pair')          
end

%% default parameters

% linewidths for plotting rectangles
plot_params.f1 = 9;    % black (background)
plot_params.f2 = 4;    % yellow (foreground)

im_path = '../Data/images/';
pt_path = '../Data/ObservedImageCoordinates/';

% control parameters for structure tensor
structure_tensor_params.tau = 0.7;   % differentiation scale
structure_tensor_params.sigma = 3;   % integration scale for structure tensor
structure_tensor_params.integration_width = ceil(1.5*structure_tensor_params.sigma)*2+1;
structure_tensor_params.sigma_n = 0.03; % noise standard deviation

%% initiate SUGR
sugr_INIT

%% read images
Image_l = imread(fullfile(im_path,Image_name_l));
Image_r = imread(fullfile(im_path,Image_name_r));
coordinate_name = fullfile(pt_path,Image_name_l(1:end-4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% show images
f1 = figure('name','Left image');
imshow(Image_l);
set(f1,'Position',[0 ss(2)/2.4 ss(1)/2 ss(2)/2]);

f2 = figure('name','Right image');
imshow(Image_r);
set(f2,'Position',[ss(1)/2 ss(2)/2.4 ss(1)/2 ss(2)/2]);

if readX == 0   % measure image correspondencies
    
    [x1, x2] = measure_homologeous_Points(f1, f2, Image_l, Image_r, 8);
    
    if x1 == -1
        error('Not enough points. We need at least 8 homologeous points!')
        return                                                             %#ok<UNRCH>
    end

    X = [x1', x2'];
    
    save(coordinate_name,'X');
    
else  % read measured image correspondencies from file
    
    load(coordinate_name,'X');
    
    N = size(X,1);                                                          %#ok<NODEF>
    
    figure(f1)
    hold on    
    for n = 1:N
        plot_square_with_background(X(n,1),X(n,2),100,plot_params.f1,plot_params.f2);
    end
    
    figure(f2)
    hold on    
    for n = 1:N
        plot_square_with_background(X(n,3),X(n,4),100,plot_params.f1,plot_params.f2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate F

% Stochastic modell of observed image points
Cxx = [ ones(8,1),zeros(8,4),ones(8,1),...
        zeros(8,4),ones(8,1),zeros(8,4),ones(8,1) ]...
        /12;

% estimate F from measured images coordinates    
[F,Cff,el,Cll,er,Crr] = sugr_estimate_F_and_epipoles_algebraically(X,Cxx);
%[F,Cff,el,er] = sugr_estimate_F_and_epipoles_algebraically_n(X,Cxx);

%% start showing epipolar lines in Figure 4 for point in Figure 3

%% plot epipole in left image
f3 = figure('name','Left image');
imshow(Image_l);
set(f3,'Position',[0 0 ss(1)/2 ss(2)/2]);
hold on
plot_epipole(f3,el,Cll,plot_params)
%% plot epipole in right image
f4 = figure('name','Right image');
imshow(Image_r);
set(f4,'Position',[ss(1)/2 0 ss(1)/2 ss(2)/2]);
hold on
plot_epipole(f4,er,Crr,plot_params)


while true 
    %% wait until a button is pressed
    disp('Measure point in Figure 3 to see epipolar line in Figure 4.')
    disp('Press right mouse button to exit')
    

    %% measure poin in left image
    x = measure_Point(f3, Image_l);
    
    if x == -1
        title('Interactively stopped')
        disp('Interactively stopped')
        break
    end
    
    [im_window, success] = get_image_window(Image_l, x, structure_tensor_params.integration_width);
    
    if success==1
        
        if cov_im == 0
            Cov = eye(2)/12;
            std_xy = sqrt(diag(Cov)');
        else
            T   = structure_tensor(im_window,...
                structure_tensor_params.tau,...
                structure_tensor_params.sigma);
            
            Cov = structure_tensor_params.sigma_n^2 * ...
                inv( T * structure_tensor_params.integration_width^2 ...
                + eye(2)*10^(-6) );                                        %#ok<MINV>
            std_xy = sqrt(diag(Cov)');
        end
        disp(['standard  deviations xy in pixel: ', num2str(std_xy)])        
        
        %plot epipole into left image again
        plot_epipole(f3,el,Cll,plot_params)        
                
        %% plot ellipse in left image
        plot_ellipse_mc(x, Cov*plot_params.magnification_m^2,'-k',plot_params.f1);
        plot_ellipse_mc(x, Cov*plot_params.magnification_m^2,'-y',plot_params.f2);
        %hold off

        %% refresh image right
        figure(f4)
        imshow(Image_r);
        hold on
        
        
        % epipolar line
        lss =  F'*[x;1];
        
        % covariance matrix of lss due to uncertainty of F and point
        J = kron(eye(3),[x;1]');
        Cllh = J*Cff*J' + F'*[Cov,zeros(2,1);zeros(1,3)]*F;
        
        % plot epipolar line into right image
        epi = sugr_Line_2D(lss,Cllh);
        sugr_plot_Line_2D(epi,'-k','-k',plot_params.f1,plot_params.magnification);
        sugr_plot_Line_2D(epi,'-y','-y',plot_params.f2,plot_params.magnification);
        
        % plot epipole into right image again
        plot_epipole(f4,er,Crr,plot_params)
        %sugr_plot_Point_2D(er,'.k','-m',3,plot_params.magnification);
        
        hold off     

    end
end

return


function plot_epipole(f,e,Cee,plot_params)

figure(f)

if abs(e(3)) > 0.0001
    
    % Plot standard ellipse of epipole
    ers=sugr_Point_2D(e,Cee);
    
    sugr_plot_Point_2D(ers,'.k','-w',plot_params.f1,plot_params.magnification)
    sugr_plot_Point_2D(ers,'.k','-m',plot_params.f2,plot_params.magnification)
    
    % plot text
    [e_E,C_EE]=sugr_get_Euclidean_Point_2D(ers);   
    text(e_E(1)+ 2*sqrt(C_EE(1,1))*plot_params.magnification,...
         e_E(2)+ 2*sqrt(C_EE(2,2))*plot_params.magnification,...
        'Epipole','FontSize',14,'FontWeight','bold','Color','m')
    
end
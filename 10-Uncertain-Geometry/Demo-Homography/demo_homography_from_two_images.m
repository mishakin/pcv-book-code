%% demo_function_homography_from_two_images: 
% determination of homography between two image pairs
%
% see PCV Sect. 10.3, Fig. 10.18
%
% demo_homography(Image_name_l,Image_name_r,magnification,
%               readX,
%               grid_out,
%               point_transfer)
%
% Image_name_l,Image_name_r = input images (left_image.ext,right_imga.ext)
%                             imgages may be the same
% magnification = magnification factor for standard ellipses
% readX  = 0 Control points are determined, interactively, 
%           Points need to be given in the following sequence:
%           for each point: left image (figure1), right image (figure 2)
%           sequence = around a quadrangle (not following a diagonal), e.g.
%           1 ---- 2
%           |      |
%           4------3
%       ~= 0 read in of control points, read from file left_imag.mat
% grid out = [Nx,xa,xe,Ny,ya,y] w.r.t. unit square 00 01 11 10
%            e.g. [3,0,1,3,0,1]   -> 3x3 grid covering unit square [0,1]^2
%            e.g. [4,-1,2,5,-1,2] -> 4x4 grid covering square [-1,2]^2
% point_transfer = 0 not demo of point transfer
%                ~=0 demo of point transfer:
%                    provide point in left image
%                    transferred point shown in right image with confidence
%                    ellipse
%
% stop interaction with pressing key, or Control C on console
%
% Samples with prepared homographies; examples = 1,2
%
% Wolfgang Förstner 3/2011
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de
%
% ........................................................................

clearvars
close all

% choose example
% example = 1; % Busstation
example = 2; % Facade

addpath(genpath('../../General-Functions'))
addpath(genpath('../Functions'))

global print_option_estimation

%% Initialize
sugr_INIT;
print_option_estimation = 0;

switch example
    case 1
        Image_name_l = 'Images/IMG_1067-Hamburg-04.JPG';
        Image_name_r = 'Images/IMG_1065-Hamburg-04.JPG';
        magnification = 3;
        readX = 1;
        grid_out = [3,0,1,3,0,1];
        point_transfer = 1;
    case 2
        Image_name_l = 'Images/IMG_1209-sect-1.JPG';
        Image_name_r = 'Images/IMG_1209-sect-1.JPG';
        magnification = 3;
        readX = 1;
        grid_out = [4,-1,2,4,-1,2];
        point_transfer = 1;
end

[Image_l,Image_r] = h_init(Image_name_l,Image_name_r);

%% determine H
[H,X] = determine_Homography_2D_from_two_images(Image_l,Image_r,readX);
close all

%% Transfer grid 
h_transfer_grid(H,X,Image_l,Image_r,magnification,grid_out);

%% Transfer single points
if point_transfer ~= 0
    h_transfer_points(H,X,Image_l,Image_r,magnification);
end





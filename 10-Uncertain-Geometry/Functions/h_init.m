% demo_h_init: initialize: initialize parameters and names for demo_homography
%
% Wolfgang Förstner 3/2011
% wfoerstn@uni-bonn.de

function [Image_l,Image_r] = h_init(Image_name_l,Image_name_r)

global tau
global sigma
global width
global sigma_n
global coordinate_sigma_x
global coordinate_name

close all

% initialize parameters for structure tensor
tau = 0.7;
sigma = 3;
width = 3*sigma;
sigma_n = 1;

maxI = 640;

Image_l = imread(Image_name_l);
Image_r = imread(Image_name_r);

[Nl,Ml,~] = size(Image_l);
factor = min(maxI/Nl,maxI/Ml);
Image_l = imresize(Image_l,factor);

[Nr,Mr,~] = size(Image_r);
factor = min(maxI/Nr,maxI/Mr);
Image_r = imresize(Image_r,factor);

coordinate_name = Image_name_l(1:end-4);
coordinate_sigma_x = 0.5;

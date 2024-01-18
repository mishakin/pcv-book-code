%% sugr_INIT
% 
% set global variables
%
% wf 1/2011
% sw 6/16 added global param for numerical precision

% function sugr_INIT

% names of types of geometric elements and transformations
global type_name

% numericla precision for comparing elements
global Threshold_numerical_precision

% for checking homogeneous part
global Threshold_Euclidean_Normalization

% radus of disk around origin for qzuasi affine transformations 
global Threshold_preimage_line_infinity


% minimal redundancy for estimating sigma_0
global min_redundancy

% option for printing and plotting
% 0=non, 1=minimum, 2=average, 3=maximum (including tests)
%results of estimations
global print_option_estimation
% plots
global plot_option
% print in general
global print_option


%% Entities
type_name{1} = 'Point_2D';
type_name{2} = 'Line_2D';
type_name{3} = 'Point_3D';
type_name{4} = 'Line_2D';
type_name{5} = 'Plane';
type_name{6} = 'Conic';
type_name{7} = 'Quadric';
type_name{8} = 'Point_Pair_2D';
type_name{9} = 'Circle';
type_name{10} = 'Ellipse';

type_name{11} = 'Point_Pair_3D_3D';
type_name{12} = 'Point_Pair_3D_2D';
type_name{12} = 'Plane_Pair';


%% Transformations 2D
type_name{20} = 'Homography_2D';
type_name{21} = 'Affinity 2D';
type_name{22} = 'Similarity 2D';
type_name{23} = 'Motion 2D';

type_name{25} = 'Essential_Matrix';
type_name{27} = 'Fundamental_Matrix';

%% Transformations 3D
type_name{30} = 'Homography_3D';
type_name{31} = 'Affinity 2D';
type_name{32} = 'Similarity 2D';
type_name{33} = 'Motion 2D';

%% Projections
type_name{123} = 'Point_Projection_2D_1D';
type_name{134} = 'Point_Projection_3D_2D';
type_name{143} = 'Point_Projection_2D_3D';
type_name{136} = 'Line_Projection_3D_2D';
type_name{163} = 'Line_Projection_2D_3D'; 

%% Korrelations
type_name{205} = 'Essential_Matrix';
type_name{207} = 'Fundamental_Matrix';

Threshold_Euclidean_Normalization = 10^(-8);
Threshold_preimage_line_infinity  = 100;
Threshold_numerical_precision = 10^(-12);
print_option_estimation = 2;
min_redundancy = 30;

plot_option = 2;
print_option = 2;


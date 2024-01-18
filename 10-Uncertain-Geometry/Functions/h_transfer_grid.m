% demo_h_transfer_grid: transfer point in demo_homography
%
% demo_h_transfer_grid(H,X,Image_l,Image_r,magnification,grid_inp,grid_out);
%
% grid_inp = 4 x 2 array: x',y' coordinates in left images 00 01 11 10
% grid_out = .Nx,.xa,.xe,.Ny,.ya,.ye of unit square 00 01 11 10
%
% Wolfgang Förstner 7/2011
% wfoerstn@uni-bonn.de

function h_transfer_grid(H,CP,Image_l,Image_r,magnification,grid_out)

ss = plot_init;

global coordinate_sigma_x

% left image
figure('Color','w','Name', 'left image','Position',[50,100,0.45*ss(1),0.45*ss(1)]);
imshow(127+0.5*Image_l);
hold on
xy = CP;

% reference quadrangle 
for n=1:4
    plot_square_with_background(xy(n,1),xy(n,2),10);
    plot_square_with_background(xy(n,1),xy(n,2),10);
end

plot([xy(1,1),xy(2,1)],[xy(1,2),xy(2,2)],'-y','LineWidth',3);
plot([xy(2,1),xy(3,1)],[xy(2,2),xy(3,2)],'-y','LineWidth',3);
plot([xy(3,1),xy(4,1)],[xy(3,2),xy(4,2)],'-y','LineWidth',3);
plot([xy(4,1),xy(1,1)],[xy(4,2),xy(1,2)],'-y','LineWidth',3);


plot([xy(1,1),xy(2,1)],[xy(1,2),xy(2,2)],'-k','LineWidth',1);
plot([xy(2,1),xy(3,1)],[xy(2,2),xy(3,2)],'-k','LineWidth',1);
plot([xy(3,1),xy(4,1)],[xy(3,2),xy(4,2)],'-k','LineWidth',1);
plot([xy(4,1),xy(1,1)],[xy(4,2),xy(1,2)],'-k','LineWidth',1);

% homography from unit square to reference quadrangle
P.h = [...
    0 0 1 xy(1,1) xy(1,2) 1;...
    0 1 1 xy(2,1) xy(2,2) 1;...
    1 1 1 xy(3,1) xy(3,2) 1;...
    1 0 1 xy(4,1) xy(4,2) 1;...
    ];
P.Crr = zeros(4,4,4);
H1 = sugr_estimation_algebraic_Homography_2D_from_point_pairs(P);


% define grid
Nx = grid_out(1);
xa = grid_out(2);
xe = grid_out(3);
dx = (xe-xa)/(Nx-1);
Ny = grid_out(4);
ya = grid_out(5);
ye = grid_out(6);
dy = (ye-ya)/(Ny-1);
N = Nx*Ny;
n = 0;
C0 = coordinate_sigma_x*diag([1,1]);
%C0=1*diag([1,1]);
for x = xa:dx:xe
    for y = ya:dy:ye
        n = n+1; 
        xyt = ([x,y,1]*H1.H')';
        xye = xyt(1:2)./xyt(3);
        gr1(n) = sugr_Point_2D(xye,C0);                                    %#ok<*AGROW>
        sugr_plot_Point_2D(gr1(n),'.k','-k',6,magnification);
        sugr_plot_Point_2D(gr1(n),'oy','-y',2,magnification);
    end
end

% transform grid an plot into 2. image with ellipse
figure('Color','w','Name', 'right image','Position',[100+0.45*ss(1),100,0.45*ss(1),0.45*ss(1)])
imshow(127+0.5*Image_r);
hold on
 % reference quadrangle 
d=2; 
for n=1:4
    plot_square_with_background(xy(n,1+d),xy(n,2+d),10);
    plot_square_with_background(xy(n,1+d),xy(n,2+d),10);
end

plot([xy(1,1+d),xy(2,1+d)],[xy(1,2+d),xy(2,2+d)],'-k','LineWidth',6);
plot([xy(2,1+d),xy(3,1+d)],[xy(2,2+d),xy(3,2+d)],'-k','LineWidth',6);
plot([xy(3,1+d),xy(4,1+d)],[xy(3,2+d),xy(4,2+d)],'-k','LineWidth',6);
plot([xy(4,1+d),xy(1,1+d)],[xy(4,2+d),xy(1,2+d)],'-k','LineWidth',6);


plot([xy(1,1+d),xy(2,1+d)],[xy(1,2+d),xy(2,2+d)],'-y','LineWidth',2);
plot([xy(2,1+d),xy(3,1+d)],[xy(2,2+d),xy(3,2+d)],'-y','LineWidth',2);
plot([xy(3,1+d),xy(4,1+d)],[xy(3,2+d),xy(4,2+d)],'-y','LineWidth',2);
plot([xy(4,1+d),xy(1,1+d)],[xy(4,2+d),xy(1,2+d)],'-y','LineWidth',2);


for n=1:N
    % transform with uncertainty
    gr2(n) = sugr_transform_with_Homography_2D(H,gr1(n));
    %% plot ellipse in right image
    sugr_plot_Point_2D(gr2(n),'ow','-w',4,magnification);
    sugr_plot_Point_2D(gr2(n),'.k','-k',2,magnification);
end
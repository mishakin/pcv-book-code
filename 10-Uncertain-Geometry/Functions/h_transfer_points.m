% demo_h_transfer_points: transfer point in demo_homography
%
% Wolfgang Förstner 3/2011
% wfoerstn@uni-bonn.de

function h_transfer_points(H,CP,Image_l,Image_r,magnification)

figure(1)
imshow(Image_l);
hold on
h_plot_CP(CP,1);

while true
    %% wait until a button is pressed
    figure(1)
    disp('grab point with mouse in image 1 or press keyboard for stop')
    title('grab point with mouse or press keyboard for stop')
    keydown = waitforbuttonpress;    
    % if a keyboard button was pressed, stop cycle
    
    if keydown == 1
        break
    end
    %% grab mouse position
    figure(1)
    imshow(Image_l);
    hold on
    h_plot_CP(CP,1);
    X=get(gca,'currentPoint');
    
    %% refresh image left
    % save only coordinates
    z = X(1,1:2)';
    
    %% plot ellipse in left image
    plot_square_with_background(z(1),z(2),10);
    hold on
    xl = sugr_Point_2D(z,10^(-4)*eye(2));
    
    %% refresh image right
    %
    figure(2)
    imshow(Image_r);
    hold on
    h_plot_CP(CP,2);
    % transferred point
    
    xr = sugr_transform_with_Homography_2D(H,xl);
    
    %% plot ellipse in right image
    sugr_plot_Point_2D(xr,'.y','-y',3,magnification);
    sugr_plot_Point_2D(xr,'ok','-k',1,magnification);
    hold off
end


% get_image_window: cuts image subwindow from large one at given position
%
% [window,success] = get_image_window(Image,z,width)
%
% Wolfgang Förstner 2/2011
% wfoerstn@uni-bonn.de

function [window,success] = get_image_window(Image,z,width)

[N,M,~]= size(Image);

% width = odd
hwidth = (width-1)/2;

y = round(z(1));
x = round(z(2));

if x-width > 0 && y-width > 0 && x < N-width && y < M-width
    window = Image(x-hwidth:x+hwidth, y-hwidth:y+hwidth,:);
    success = 1;
else
    success = 0;
    window  = 0;
end

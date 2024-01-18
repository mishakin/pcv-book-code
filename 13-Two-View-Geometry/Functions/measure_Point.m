function [x] = measure_Point(f, I, N, cs)
% This function enables the user to measure points in an image
%
% Inputs:
%   f - figure handle
%   I - image
%   N - number of points to measure, default 1
%   cs - coordinate system
%        'matlabimagecs' - default, Matlab image coordinatesystem
%                          (0,0) top left
%                          1.dim to the left, 2.dim top down
%        'xy'              righthandside coord.system
%                          (0,0) top left
%                          1.dim top down, 2.dim to the left
%        'xy_bl'           righthandside coord.system
%                          (0,0) bottom left
%                          1.dim to the left, 2.dim bottom up
%
% Outputs:
%   x - measured coordinates
%       dimension 2xN
%
% Author:
%   Falko Schindler (falko.schindler@uni-bonn.de)
%   Susanne Wenzel (swenzel@igg.unibonn.de)
%
% Date:
%   December 2010
% Last change
%   April 2018 Susanne Wenzel  (swenzel@igg.unibonn.de)

if nargin<4
    cs = 'matlabimagecs';
end
if nargin<3
    N = 1;
end

figure(f) 

%% measure coordinates
x = zeros(2,N);

b = 100;
for n = 1 : N   
       
    title('Identify rough position for zoom, or right klick to exit.')
    hold off
    
    [x0,y0,button] = ginput_10(1);
    if button == 3
        x = -1;
        return 
    end
    
    xrange = round(max(x0 - b, 1) : min(x0 + b, size(I, 2)));
    yrange = round(max(y0 - b, 1) : min(y0 + b, size(I, 1)));
    
    imshow(I(yrange, xrange, :), ...
        'Border', 'tight', 'InitialMagnification', 'fit', ...
        'XData', xrange, 'YData', yrange);
    title('Measure image point')
    
    xy = ginput_10(1);
    
    switch cs
        case 'matlabimagecs'
            x(:, n) = [xy(1); xy(2)];
        case 'xy'
            x(:, n) = [xy(2); xy(1)];
        case 'xy_bl'
            x(:, n) = [xy(1); size(I, 1)-xy(2)+1];
        otherwise
            error('Unknown coordinate system')            
    end
        
end

imshow(I, 'Border', 'tight', 'InitialMagnification', 'fit');
hold on

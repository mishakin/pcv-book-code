function [x1, x2] = measure_homologeous_Points(f1, f2, I1, I2, N, cs)
% This function enables the user to measure homologous points in two
% images.
%
% Inputs:
%   f1, f2 - figure handles
%   I1, I2 - images
%   N - number of points to measure
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
%   - x1 and x2  measured coordinates
%       dimension 2xN each
%
% Author:
%   Falko Schindler (falko.schindler@uni-bonn.de)
%   Susanne Wenzel (swenzel@igg.unibonn.de)
%
% Date:
%   December 2010
% Last change
%   April 2018 Susanne Wenzel  (swenzel@igg.unibonn.de)


if nargin<6
    cs = 'matlabimagecs';
end
if nargin<5
    N = 1;
end

%% measure coordinates
x = {zeros(2, N);
     zeros(2, N)};
 
I = {I1;I2};
f = {f1;f2};

b = 100;
for n = 1 : N
    for i = 1 : numel(I)
        if n>1 && i==1
            close(f_tmp{:})
        end
        figure(f{i})  
        hold off     
        title('Identify rough position for zoom, or klick right to exit.')
        disp(['Identify rough position for zoom in image ',num2str(i),' or right klick to exit.'])
        
        [x0,y0,button] = ginput_10(1);
        if button == 3
            x1 = -1;
            x2 = -1;
            return 
        end
        
        xrange = round( max(x0 - b, 1) : min(x0 + b, size(I{i}, 2)) );
        yrange = round( max(y0 - b, 1) : min(y0 + b, size(I{i}, 1)) );
        
        imshow(I{i}(yrange, xrange, :), ...
            'Border', 'tight', 'InitialMagnification', 'fit', ...
            'XData', xrange, 'YData', yrange);
        
        title('Measure exact position.')
        disp(['Measure exact position in image ',num2str(i)])
        xy = ginput_10(1);        
        
        switch cs
            case 'matlabimagecs'
                x{i}(:, n) = [xy(1); xy(2)];
            case 'xy'
                x{i}(:, n) = [xy(2); xy(1)];
            case 'xy_bl'
                x{i}(:, n) = [xy(1); size(I, 1)-xy(2)+1];
            otherwise
                error('Unknown coordinate system')            
        end
        
        imshow(I{i}, 'Border', 'tight', 'InitialMagnification', 'fit');   
        hold on
        for m = 1:n
            plot_square_with_background(x{i}(1, m),x{i}(2, m),50,8,4);
        end
        
        fp = get(gcf,'Position');
        ss = get(0,'ScreenSize');
        if ss(4)-(fp(4)+fp(2)) < fp(2)
            % show zoom window below current figure
            fs = ss(4)-fp(4);
            f_tmp{i} = figure('name','Last point, detail','Position',[(fp(1)+fp(3))/2-fs/2, 10, fs, 2/3*fs]); %#ok<*AGROW>
        else
            fs = ss(4)-(fp(4)+fp(2));
            f_tmp{i} = figure('name','Last point, detail','Position',[(fp(1)+fp(3))/2-fs/2, ss(4)-fs, fs, fs]);
        end
        imshow(I{i}(yrange, xrange, :))
        hold on
        plot_square_with_background(xy(1)-min(xrange), xy(2)-min(yrange),50,8,4);
    end
    
    figure(f{1})  
    disp('For next point click into image 1')
    w = waitforbuttonpress;
    if w == 0
        continue
    else
        continue
    end
end
    
[x1, x2] = deal(x{:});



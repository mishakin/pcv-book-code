%% plot_init: initialize visualization
%
% usage: plot_init
%             set Latex as default interpreter for all text props
%        screensize = plot_init   
%             set Latex as default interpreter for all text props and
%             return screensize: 1x2 int, size of current screen [width, height]
%
% Susane Wenzel 09/16
% wenzel@igg.uni-bonn.de

function varargout = plot_init

if verLessThan('matlab', '8')
else
%     set Latex as default interpreter for all text props
    set(0, 'defaulttextinterpreter','latex');
    set(0, 'defaultAxesTickLabelInterpreter','latex');
    set(0, 'defaultLegendInterpreter','latex');
end

if nargout==1
    % get current screensize, for proper positioning of figures
    screensize = get(0,'ScreenSize'); screensize = screensize(3:4);
    varargout = {screensize};
end


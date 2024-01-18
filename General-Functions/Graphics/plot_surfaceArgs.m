%% parses arguments for plotting surfaces
%
% Usage:
%    args = plot_surfaceArgs(varargin)
%
%    varargin - key-value pairs according to matlabs plot parameters
%    args - key-value pairs saved as struct
%
% Susanne Wenzel 09/16
% wenzel@igg.uni-bonn.de
function args = plot_surfaceArgs(varargin)

% set default parameters
args.axes = gca;
args.FaceLighting = 'flat';    % standard Matlab prop
args.EdgeColor = 'k';          % standard Matlab prop
args.alpha = 1;
args.view = [-20,30];
args.shading = 'none';       % standard Matlab prop
args.ColorFct = 'none';      % option: 'smoothtanh' 
args.Grid = 'off';
args.colormap = gray;
args.plotfun = @surf;

args = parseArguments(args, varargin{:});




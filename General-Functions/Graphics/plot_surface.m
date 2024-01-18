%% Ploting 2.5D surfaces
%
% Usage:
%    plot_surface(ds,BB,delta_x)
%    h = plot_surface(ds,BB,delta_x, [varargin])
% 
% ds: double NxM, grid containing surface to plot
% BB: double 1x4 Bounding Box of surface footprint
%     [xmin, ymin, xmax, ymax]
% delta_x: double, grid size
% varagin: key-value pairs according to matlabs plot parameters
%

function h = plot_surface(ds,BB,delta_x,varargin)

args = plot_surfaceArgs(varargin{:});

Nr = size(ds,1); Mc = size(ds,2);

X = ([0:Nr-1]'*ones(1,Mc))*delta_x+BB(1);                                   %#ok<*NBRAK>
Y = (ones(Nr,1)*[0:Mc-1])*delta_x+BB(2);

colormap(args.colormap)

switch args.ColorFct
    case 'smoothtanh'    
        shade = -[ 0 -1; 1 0]/sqrt(2);
        col=conv2(ds,shade,'same');
        %col=(1+col./sqrt(1+col.^2))/2;
        colo=min(1, max(0,tanh(100*(col-0)/(max(col(:))-min(col(:))))));

        h = args.plotfun(X,Y,ds,colo);
    case 'none'
        h = args.plotfun(X,Y,ds);
end
set(h,'FaceLighting',args.FaceLighting,'EdgeColor',args.EdgeColor)
set(gca,'XGrid',args.Grid)
set(gca,'YGrid',args.Grid)
set(gca,'ZGrid',args.Grid)
    
alpha(args.alpha);
if ~strcmp(args.shading,'none')
    shading(args.shading);
end
view(args.view);

%% explicit plot of conic (ellipse or hyperbola)
%
% h = plot_conic_explicit( C, ... ) plots the conic defined by the 3x3 matrix
% C explicitly. 
%
% * C = conic
% * center_type = string specifying plot type of center
% * bound_type = string specifying plot type of conic
% * linewidth = width of line of conic
% * factor = magification of conic compared to position
% * se = [start,end] of hyperbola for line segments
%        se=[0,0]'  then line segment (+-sqrt(2))
%        se=[1,1]'  then line until xlim,ylim (default)
%
% $Log: plot_conic_explicit.m,v $
% Revision 1.1  2009/08/03 11:10:33  Meidow
% *** empty log message ***
%
% adapted
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de

function sugr_plot_conic_explicit...
    (C,center_type,bound_type,linewidth,factor,se)

warning off

if nargin < 4
    linewidth=2;
end
if nargin < 5
    factor=1;
end
if nargin < 6
    se = [1,1];
end

C = C/norm(C+eps);
N0=50;
N = 2*N0;  % number of supporting points


Chh = C(1:2,1:2);
ch0 = C(1:2,  3);
c00 = C(  3,  3);

d = det( Chh +eps*eye(2));
x0 = -inv(Chh+eps*eye(2)) *ch0;            % centre point
plot( x0(1), x0(2), center_type);

c00q = c00 -ch0'*inv(Chh+eps*eye(2))*ch0;                                  %#ok<MINV>
[Rho,Lambda] = eig( -Chh/ (c00q+eps) );

lambda = diag(Lambda);

hold on
switch sign(d)
    case +1
        %% Ellipse
        plot( x0(1), x0(2), center_type);
        %  (1) plot ellipse
        if lambda(1)<0
            lambda = -lambda;
        end
        
        lambda=lambda/factor^2;

        t = linspace(0,2*pi, N);

        x = sqrt(1/lambda(1)) * sin(t);
        y = sqrt(1/lambda(2)) * cos(t);
        x = bsxfun(@plus, Rho * [x; y], x0);
        plot( x(1,:), x(2,:), bound_type,'LineWidth',linewidth);

    case -1
        %% Hyperbola
        % (2) plot hyperbola
        if lambda(2)<0
            lambda = flipud(lambda);
            Rho = fliplr(Rho);
        end
                
        lambda(2)=lambda(2)/factor^2;

        %         t = linspace(-asinh(a), +asinh(a), N);
        %         % t = linspace(-pi,pi, N);
        %
        %         x = sqrt(1/-lambda(1)) *sinh(t);
        %         y = sqrt(1/+lambda(2))* cosh(t);       
         
        
        if se(1) ~= 0 && se(2) ~= 0 || se(1) ~= 1 && se(2) ~= 1     % bounds given
            t = linspace(atan(se(1)*sqrt(-lambda(1))),atan(se(2)*sqrt(-lambda(1))),N);
        end
        if se(1) == 1 && se(2) == 1            % bounds by xlim and ylim
            dir = Rho(:,1);
            xl  = xlim;
            yl  = ylim;
            s1  = dir'*[xl(1)-x0(1);yl(1)-x0(2)];
            s2  = dir'*[xl(1)-x0(1);yl(2)-x0(2)];
            s3  = dir'*[xl(2)-x0(1);yl(2)-x0(2)];
            s4  = dir'*[xl(2)-x0(1);yl(1)-x0(2)];
            se(1) = min([s1,s2,s3,s4])-1;
            se(2) = max([s1,s2,s3,s4])+1;
            t = linspace(atan(se(1)*sqrt(-lambda(1))),atan(se(2)*sqrt(-lambda(1))),N);
        end
        if  se(1) == 0 && se(2) == 0                      % line segment 
            t = linspace(-pi/4,pi/4,N);
        end
        x = sqrt(1/-lambda(1)) * tan(t);
        y = sqrt(1/+lambda(2)) ./ cos(t);

        x1 = repmat(x0, 1, N) + Rho*[x; +y];
        x2 = repmat(x0, 1, N) + Rho*[x; -y];
        plot( x1(1,:), x1(2,:), bound_type,'LineWidth',linewidth);
        plot( x2(1,:), x2(2,:), bound_type,'LineWidth',linewidth);
        plot( [(x1(1,1)+x2(1,1))/2,(x1(1,N)+x2(1,N))/2],...
            [(x1(2,1)+x2(2,1))/2,(x1(2,N)+x2(2,N))/2] ,center_type);


    otherwise
        error('degenrated conic')
end

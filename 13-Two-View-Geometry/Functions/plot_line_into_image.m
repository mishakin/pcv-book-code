% plot a line into an image if possible with error bound
%
% Image  = Image
% line   = homogeneous line coordinates
% C      = covariance matrix of homogeneous coordinates
%
% plot_line_into_image(Image_r,lss);
%
% Wolfgang Förstner 10/2010
% wfoerstn@uni-bonn.de

function plot_line_into_image(Image,line,C,f1,f2)

% default widths
if nargin < 4
    f1 = 4;
    f2 = 2;
end

% normalized normal vector
n = line(1:2)/norm(line(1:2));

% size of image
[N,M,~] = size(Image);

% plot clipped line
[x1,y1,x2,y2] = clipline(line(1:2),line(3),[0,M,0,N]);
plot([x1,x2],[y1,y2],'-k','LineWidth',f1);
plot([x1,x2],[y1,y2],'-y','LineWidth',f2);

if C(1,1) ~= 0
    % plot hyperbolic band
    L = 30;
    alpha = 1:-1/(L-1):0;

    % initiate reference point
    xah = [x1,y1,1]';
    sqa = sqrt(xah'*C*xah/norm(line(1:2))^2);
    % left and right point
    xal = [x1;y1]+n*sqa;
    xar = [x1;y1]-n*sqa;
    for l = 1:L-1
        % 2. reference point
        xe = alpha(l+1)*x1+(1-alpha(l+1))*x2;
        ye = alpha(l+1)*y1+(1-alpha(l+1))*y2;
        xeh = [xe,ye,1]';
        sqe = sqrt(xeh'*C*xeh/norm(line(1:2))^2);
        % left right point
        xel = [xe;ye]+n*sqe;
        xer = [xe;ye]-n*sqe;
        % left right line

        % plot clipped line left and right
        plot([xal(1),xel(1)],[xal(2),xel(2)],'-k','LineWidth',f1);
        plot([xal(1),xel(1)],[xal(2),xel(2)],'-y','LineWidth',f2);
        plot([xar(1),xer(1)],[xar(2),xer(2)],'-k','LineWidth',f1);
        plot([xar(1),xer(1)],[xar(2),xer(2)],'-y','LineWidth',f2);

        % transfer
        % left and right point
        xal = xel;
        xar = xer;

    end
end

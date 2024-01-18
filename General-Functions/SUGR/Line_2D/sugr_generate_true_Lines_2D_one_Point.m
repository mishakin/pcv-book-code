% sugr_generate_true_2d_lines_one_point: generates image lines through one point
%
% lines = sugr_generate_true_2d_lines_one_point(x0,d,N,lengthq,resolution)
%
% * x0 true intersection point, homogeneous vector, direction from projection centre
% * d = image size d*[-1,1] x  d*[-1,1] ---- assuming c=1
% * N  = total number of lines to be generated
% * lengthq [pel] = average length of lines in pixels
%
% *     lines have length lengthq*f, log(f) in [-1,3],
% *     thus f in[1/e, e^3] appr. [0.3,25]
%
% * resolution [pel] = number of pixels (square image)
% 
% lines = N x 4 array with [u0,v0,phi_0,lengthl]    
%
% Wolfgang Förstner 2/2011
% wfoerstn@uni-bonn.de
% 
% see also sugr_Line_2D, sugr_perturb_Lines_2D

function lines = sugr_generate_true_Lines_2D_one_Point(x0,d,N,lengthq,resolution)


%% Initiate
%
% determine given point

c=1; % nur zum Test auf Sensitivität gegen Änderungen von c

van_1 = x0/norm(x0);
van_1(3) = van_1(3)*c;

%% prepare generation of lines

lines        = zeros(N,4);

ii = 0;

% Delta_x=d/resolution;

maxl = 0;
minl = 10000;

%% generate N lines
for n = 1:N
    % generate line to point van_1

    % Generate true data
    
    % centre, randomly distributed in image
    z0     = (rand(1)*2-1+1i*(rand(1)*2-1))*d;

    % true line
    e_line = cross([real(z0),imag(z0),1],van_1);
    if rand(1) < 0.5 e_line=-e_line; end;

    % true direction of line
    phi_0 = atan2(e_line(2),e_line(1))+pi/2;

    % length
    log_f=(rand(1)*4-1)/2;
    f=exp(log_f);
    lengthl= lengthq*f;
    minl=min(minl,lengthl);
    maxl=max(maxl,lengthl);
    ii=ii+1;
    lines(ii,:)=[real(z0),imag(z0),phi_0,lengthl];
end;

% minlength_maxlength_of_segments=[minl,maxl];
disp(['Minimum and maximum length of segments : [', num2str(minl),',',num2str(maxl),']']);


            
  





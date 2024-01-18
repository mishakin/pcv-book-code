%% determine structure tensor from image patch
% 
% im    = window
% tau   = differentiation scale
% sigma = integration scale
%
% T     = structure tensor 
%       = \overline(g'_x g_x)
%       = G(\sigma) * (G(\tau)*\nabla_g'  G(\tau)*\nabla_g )
%
% T     = structure_tensor(im,tau,sigma)
%
% Observe: 
% The structure tensor is not the sum ...
%                         but the average of the squared gradient
%
% Wolfgang Förstner 4/2018
% wfoerstn@uni-bonn.de

function T  = structure_tensor(im,tau,sigma)

% size of image, number of channels
[N,M,K] = size(im);

% initiate T
T = zeros(2);

% sum up over all channels
for k = 1:K
    % gradient at scale tau
    im_x = gaussFFT(double(im(:,:,k))/256,tau,'Gx'); 
    im_y = gaussFFT(double(im(:,:,k))/256,tau,'Gy'); 
    % squared gradient
    imxx = im_x.*im_x;
    imyy = im_y.*im_y;
    imxy = im_x.*im_y;
    % mean squared gradient
    txx = gaussFFT(imxx,sigma,'G');
    tyy = gaussFFT(imyy,sigma,'G');
    txy = gaussFFT(imxy,sigma,'G');
    % sum up channels
    T = T+[txx((N+1)/2,(M+1)/2) , txy((N+1)/2,(M+1)/2); ...
           txy((N+1)/2,(M+1)/2) , tyy((N+1)/2,(M+1)/2)];
end

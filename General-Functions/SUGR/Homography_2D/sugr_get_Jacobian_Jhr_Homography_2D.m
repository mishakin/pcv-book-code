%% Get Jacobian Jhr dh/dr (h = vec H) for estimating 2D Homography  
% 
% [Jhr,Jkr,Jhk] = sugr_get_Jacobian_Jhr_Homography_2D(H)
% 
% H   = 3x3-matrix, spectrally normalized Homography
%
% Jhr = 9x8-matrix Jacobian dh/dr,r -> h
% Jkr = 9x9-matrix Jacobian dk/dr,r -> k
% Jhk = 9x9-matrix Jacobian dh/dk,k -> h
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 2/2011

function [Jhr,Jkr,Jhk] = sugr_get_Jacobian_Jhr_Homography_2D(H)

Ht  = H';      
Jkr = [eye(8); [-1,0,0,0,-1,0,0,0]];    % 9x8-Jacobian dk/dr,r  -> k   
Jhk = kron(Ht,eye(3));                  % 9x9-Jacobian dh/dk,k  -> h
Jhr = Jhk * Jkr;                        % 9x8-Jacobian dh/dr,r  -> h

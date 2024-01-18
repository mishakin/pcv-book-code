%% Generate 2D point pairs from homography,
% i.e., for given H generate N point pairs
%
% PP = sugr_generate_true_2D_point_pairs_homography(H,N,boolean_r);
%
% H = 3x3 matrix
% N = number of points (in square [-1,1]^2
% boolean_r = boolean: points should sit random
%                else: points sit in a square (N should be square)
%
% PP = point pairs
%      PP.h = N x 6 matrix of pairs of homogeneous point coordinates
%      PP.Crr = N x 4 x 4 sith 4 x4 CovM of point pairs
%      PP.type = 8 * ones(N,1)
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 02/2011

function PP = sugr_generate_true_2D_point_pairs_homography(H,N,br)

if br
    for n=1:N
        x             = [rand(2,1)*2-1;1];
        y             = H * x;
        PP.h(n,:)     = [x',y'];
        PP.Crr(n,:,:) = zeros(4);
        PP.type(n)    = 8;
    end
else
    M = ceil(sqrt(N));
    n=0;
    for m=1:M
        for k=1:M
            n = n+1;
            if n <= N
                x = [-(M-1)/2+2*(m-1)/(M-1);-(M-1)/2+2*(k-1)/(M-1);1];
                y             = H * x;
                PP.h(n,:)     = [x',y'];
                PP.Crr(n,:,:) = zeros(4);
                PP.type(n)    = 8;
            end
        end
    end
end

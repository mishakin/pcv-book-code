%% create 2D homography
%
% H = sugr_Homography_2D(a1,a2);
%
% input cases
% (1) 8-vector of Euclideanly normalized homography, CovM == 0
% (1) list of Euclidean points, Nx2 matrix of x,y-coordinates -> H = conditioning matrix
% (1) 3 x 3-matrix, CovM == 0
% (1) list of structs of points and lines -> generate conditioning matrix
% (2) 3 x 3-matrix with full covariance matrix
% (2) 8-vector Euclidean homography with full covariance matrix
%
% Homography = structure
% *             .H   = 3 x 3-Matrix, determinant normalized (det= +1)
% *             .Crr = covariance matrix of 8 parameters of differential left factor
% *             .type = 20
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 1/2011

function Homography_2D = sugr_Homography_2D(a1,a2)

% default
Chh  = zeros(9);
%
switch nargin
    case 1
        % 9-matrix or 8-vector
        switch isstruct(a1)
            case 0
                [~,columns_a1] = size(a1);
                switch columns_a1
                    case 1 % 8-vector of elements except H(3,3)
                        %% (1) 8-vector of Euclidean normalized homography, CovM == 0
                        H = reshape([a1;1],3,3);
                        Homography_2D = sugr_minimal_Homography_2D(H,Chh);
                    case 2 % list of Euclidean points
                        xy = a1;
                        max_xy=max(xy);
                        min_xy=min(xy);
                        
                        H = [1 ,0 ,-(min_xy(1)+max_xy(1))/2;...
                            0 ,1 ,-(min_xy(2)+max_xy(2))/2;...
                            0 ,0 ,max(max_xy(1)-min_xy(1),max_xy(2)-min_xy(2))];
                        Homography_2D = sugr_Homography_2D(H);
                    case 3
                        %% (1) 3 x 3-matrix, CovM == 0
                        Homography_2D = sugr_minimal_Homography_2D(a1,Chh);
                        
                        
                end
            case 1
                %% (1) generate conditioning matrix
                % set of point structs and line structs -> H = conditioning matrix
                X =a1;
                N  = length(X);
                xy = zeros(N,2);
                % for list of structs
                for i=1:N
                    if X(i).type == 1
                        [xe,~] = sugr_get_Euclidean_Point_2D(X(i));
                        xy(i,:) = xe';
                    end
                    
                    if X(i).type == 2
                        f = sugr_get_centroid_Line_2D(X(i));
                        [fe,~] = sugr_get_Euclidean_Point_2D(f);
                        xy(i,:) = fe';
                    end
                end
                
                max_xy=max(xy);
                min_xy=min(xy);
                
                H = [1 ,0 ,-(min_xy(1)+max_xy(1))/2;...
                    0 ,1 ,-(min_xy(2)+max_xy(2))/2;...
                    0 ,0 ,max(max_xy(1)-min_xy(1),max_xy(2)-min_xy(2))];
                
                Homography_2D = sugr_Homography_2D(H);
        end
    case 2 %
        %% 9-matrix or 8-vector with full covariance matrix
        size_a1 = size(a1,1);
        switch size_a1
            case 3
                %% (2) 3 x 3-matrix
                Homography_2D = sugr_minimal_Homography_2D(a1,a2);
                
            case 8
                %% (2) 8-vector Euclidean homography
                H = reshape([a1;1],3,3);
                Chh = [a2 zeros(8,1);zeros(1,9)];
                Homography_2D = sugr_minimal_Homography_2D(H,Chh);
                
        end
        
end

Homography_2D.type = 20;

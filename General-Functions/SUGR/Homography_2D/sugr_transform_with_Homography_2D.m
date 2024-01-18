%% transform uncertain elements using a projective mapping
%
% y = sugr_transform_Homography_2D(H,x);
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 1/2011
% sw 9/2016

function y = sugr_transform_with_Homography_2D(H,x)

N = length(x);

switch H.type
    case 20
        switch N
            case 1
                switch x.type
                    case 1
                        %% transform a 2D point:                                              x' = H * x
                        Jx     = H.H;                                 % Jacobian x' -> x
                        y      = Jx * x.h;                            % x'
                        Jh     = kron(x.h', eye(3));                  % Jacobian x' -> h
                        CHhh   = sugr_get_CovM_homogeneous_Homography_2D(H);
                        Cxhh   = sugr_get_CovM_homogeneous_Vector(x);
                        Cyhh   = Jx * Cxhh * Jx' + Jh * CHhh * Jh';   % C_x'x_
                        y      = sugr_Point_2D(y,Cyhh);
                        
                    case 2
                        %% transform a 2D line:                                                 l' = inv(H') * l
                        Jx     = inv(H.H');                           % Jacobian l' -> l
                        ls     = Jx * x.h;                            %#ok<*MINV> % l'
                        Jh     = kron(x.h', eye(3));                  % Jacobian l' -> h
                        CHhh   = sugr_get_CovM_homogeneous_Homography_2D(H);
                        Cxhh   = sugr_get_CovM_homogeneous_Vector(x);
                        Clhh   = Jx * Cxhh * Jx' + Jh * CHhh * Jh'; % C_l'l'
                        y      = sugr_Line_2D(ls,Clhh);
                        
                end
            otherwise
                %% Transform a list of elements
                y = zeros(N,1);
                for n=1:N
                    y(n) = sugr_transform_with_Homography_2D(H,x(n));                    % transform each element
                end
        end
end


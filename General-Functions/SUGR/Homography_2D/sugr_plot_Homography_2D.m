%% plot uncertain homography 2D
%
% sugr_plot_Homography_2D(H,center_type,bound_type,lw,factor,type)
%
% * H             = homography
% * center_type   = type of plotting the element
% * bound_type    = type of plotting the confidence region
% * lw            = linewidth
% * factor        = magnification factor
% * type          = type for showing homography
%                 = 'x' only points are uncertain
%                 = 'Hx' points and homography are uncertain
%                 = 'H' only homopgraphy is uncertain
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 1/2011

function sugr_plot_Homography_2D(H,center_type,bound_type,lw,factor,type)

if nargin < 6
    type='';
end

%% plot 2D homography by original and transfomed grid

xlim([-1.4,4.3]);
ylim([-3.8,1.4]);

Ni = 3;
Nj = 3;
k  = 0;
switch type
    case 'x'
        Cxx = 0.0004*eye(2);
        CHH = zeros(9);
        Hp = sugr_minimal_Homography_2D(H.H,CHH);
    case 'Hx'
        Cxx = 0.000004*eye(2);
        Hp = H;
    case 'H'
        Cxx = eye(2)*10^(-10);
        Hp = H;
    otherwise
        Cxx = zeros(2);
        CHH = zeros(9);
        Hp = sugr_minimal_Homography_2D(H.H,CHH);
end


for i=1:Ni
    for j=1:Nj
        k = k+1;
        X(k) = sugr_Point_2D([(i-Ni/2-1/2), (j-Nj/2-1/2)]',Cxx);
    end
end
Nk    =  Ni*Nj+1;
X(Nk) = X(1);

for k=1:Nk
    Y(k) = sugr_transform_with_Homography_2D(Hp,X(k));
    Y(k).Crr;
end
if  sugr_get_isfinite_Point_2D(Y(k))
    for k=1:Nk
        sugr_plot_Point_2D(X(k),center_type,'-r',lw,factor);
        hold on
        sugr_plot_Point_2D(Y(k),center_type,'-b',lw,factor);
    end
    for i=1:Ni
        [Xs,~] = sugr_get_Euclidean_Point_2D(X(Nj*(i-1)+1));
        [Xe,~] = sugr_get_Euclidean_Point_2D(X(Nj*i));
        plot([Xs(1),Xe(1)],[Xs(2),Xe(2)],'-r','LineWidth',2);
        [Ys,~] = sugr_get_Euclidean_Point_2D(Y(Nj*(i-1)+1));
        [Ye,~] = sugr_get_Euclidean_Point_2D(Y(Nj*i));
        plot([Ys(1),Ye(1)],[Ys(2),Ye(2)],'-r','LineWidth',2);
    end
    for j=1:Nj
        [Xs,~] = sugr_get_Euclidean_Point_2D(X(j));
        [Xe,~] = sugr_get_Euclidean_Point_2D(X(j+(i-1)*Nj));
        plot([Xs(1),Xe(1)],[Xs(2),Xe(2)],'-r','LineWidth',2);
        [Ys,~] = sugr_get_Euclidean_Point_2D(Y(j));
        [Ye,~] = sugr_get_Euclidean_Point_2D(Y(j+(i-1)*Nj));
        plot([Ys(1),Ye(1)],[Ys(2),Ye(2)],'-r','LineWidth',2);
    end
    
end

xlim([-1.4,4.3]);
ylim([-4.9,1.8]);
axis equal


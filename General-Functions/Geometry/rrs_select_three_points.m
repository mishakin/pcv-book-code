%% Select 3 points for spatial resection and possibly plot
%
% p_selected = rrs_select_three_points(X,y,Tol)
% 
% X    N-struct given scene points, possibly at infinity
% y    N-struct given image directions
% Tol  Tolerance for points to be not at infinity (eg 100)
%
% X3   3x3 matrix with rows being scene points (X,Y,Z)
% y3   3x3 matrix with rows being image directions


function p_selected = rrs_select_three_points(X,y,Tol)

global plot_option
global print_option

N = size(X.h,1);
% number of trials of triplets
K=30;
dmax = 0;
for k=1:K
    pe = randperm(N);
    % determinant of directions (should be large)
    d3 = abs(det([y.h(pe(1),:);y.h(pe(2),:);y.h(pe(3),:);]));
    % no point should be far away
    X_f_1 = abs(X.h(pe(1),4)) > 1/Tol;
    X_f_2 = abs(X.h(pe(2),4)) > 1/Tol;
    X_f_3 = abs(X.h(pe(3),4)) > 1/Tol;
    if d3 > dmax && X_f_1*X_f_2*X_f_3 == 1
        three = pe;
        dmax= d3;
        pe_selected = pe;
        X3 = [X.h(pe(1),1:3)/X.h(pe(1),4);...
              X.h(pe(2),1:3)/X.h(pe(2),4);...
              X.h(pe(3),1:3)/X.h(pe(3),4)];
        y3 = [y.h(pe(1),:);...
              y.h(pe(2),:);...
              y.h(pe(3),:)];
    end
end
if dmax < 0.00001 
    display('no adequate triple found')
    return
end

if print_option > 0
    display(['Points for approximate values : ',num2str(pe_selected(1:3))])
end
if plot_option > 0
    figure
    subplot(1,2,1)
    hold on
    plot3(0,0,0,'.b','MarkerSize',15)
    for n=1:3
        plot3(X3(n,1),X3(n,2),X3(n,3),'.r','MarkerSize',15)
    end
    axis equal
    subplot(1,2,2)
    hold on
    plot3(0,0,0,'.b','MarkerSize',15)
    for n=1:3
        plot3(y3(n,1),y3(n,2),y3(n,3),'.r','MarkerSize',15)
        plot3([0,y3(n,1)],[0,y3(n,2)],[0,y3(n,3)],'-r')
    end
    axis equal
end

% selected points
p_selected=pe(1:3);

end


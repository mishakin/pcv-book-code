% plot control points
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de

function h_plot_CP(X,fi)

figure(fi)
d = 2*(fi-1);

for n=1:4
    plot_square_with_background(X(n,1+d),X(n,2+d),10);
    plot_square_with_background(X(n,1+d),X(n,2+d),10);
end

plot([X(1,1+d),X(2,1+d)],[X(1,2+d),X(2,2+d)],'-y','LineWidth',3);
plot([X(1,1+d),X(3,1+d)],[X(1,2+d),X(3,2+d)],'-y','LineWidth',3);
plot([X(1,1+d),X(4,1+d)],[X(1,2+d),X(4,2+d)],'-y','LineWidth',3);
plot([X(2,1+d),X(3,1+d)],[X(2,2+d),X(3,2+d)],'-y','LineWidth',3);
plot([X(2,1+d),X(4,1+d)],[X(2,2+d),X(4,2+d)],'-y','LineWidth',3);
plot([X(3,1+d),X(4,1+d)],[X(3,2+d),X(4,2+d)],'-y','LineWidth',3);

plot([X(1,1+d),X(2,1+d)],[X(1,2+d),X(2,2+d)],'-k','LineWidth',1);
plot([X(1,1+d),X(3,1+d)],[X(1,2+d),X(3,2+d)],'-k','LineWidth',1);
plot([X(1,1+d),X(4,1+d)],[X(1,2+d),X(4,2+d)],'-k','LineWidth',1);
plot([X(2,1+d),X(3,1+d)],[X(2,2+d),X(3,2+d)],'-k','LineWidth',1);
plot([X(2,1+d),X(4,1+d)],[X(2,2+d),X(4,2+d)],'-k','LineWidth',1);
plot([X(3,1+d),X(4,1+d)],[X(3,2+d),X(4,2+d)],'-k','LineWidth',1);

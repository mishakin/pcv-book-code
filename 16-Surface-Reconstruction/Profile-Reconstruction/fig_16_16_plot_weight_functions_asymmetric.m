%% Fig. 16.16 page 756
% compare weight functions for a/symmetric outliers
%
% Wolfgang Förstner 08/14
% last changes: Susanne Wenzel 09/16
% wfoerstn@uni-bonn.de, wenzel@igg.uni-bonn.de

close all


addpath(genpath('../../General-Functions'));

%% plot settings
ss = plot_init;

%% set parameters

% resolution
MK = 1000;
M = 50;
k = 1;                % critical value
g = 2*k;              % g-parameter for Kraus (taken negative)
w = 2*k;              % w-parameter for Kraus

% set plot range
vK = -5:10/(MK-1):5;   % for Kraus
v =  -5:10/(M-1):5;    % others

%% determine Kraus' weight functions
weights = zeros(MK,1);
for m = 1:MK
    
    if vK(m) < g-w
        weights(m) = 0;
        
    elseif vK(m) > g
        weights(m) = 0.985;
        
    else
        weights(m) = 0.985/(1+(vK(m)-g)^4);
        
    end
end
w_kraus=weights;

%% determine asymmetric L_1 weight function
weights = max( min( 1,1./abs(v(1:M))/k+0.0001 ),...
               ( (sign(v(1:M))/k)+1 )/2 ...
           );
w_L1 = weights;

x = vK;
y = w_kraus';
rangeX = max(x) - min(x);
rangeY = max(y) - min(y);
Q = [diff(x) / rangeX; diff(y) / rangeY];
L = [0, cumsum(sqrt(sum(Q .* Q)))];
L = (length(L) / L(end)) * L;
a = 1;  % start point
b = 10; % gap width
c = 60; % gap frequency
index  = rem(round(L) + a, c) <= b;

%% plot
figure('name','Fig. 16.16 Weight functions for outliers','color','w',...
    'Position',[0.2*ss(1),0.2*ss(2),0.4*ss(1),0.4*ss(2)]);
hold on
%plot(v,w_kraus,'--r','LineWidth',4)
%dashline(v,w_kraus,10,2,10,2,'color','r','LineWidth',4)
plot(v,w_L1,'-b','LineWidth',6)
plot(v,min(0.985,0.95./abs(v(1:M))/k+0.0001)-0.003,'.k','MarkerSize',20)

yDashed = y;
yDashed(index) = NaN;
plot(x, yDashed,'color','r','LineWidth',6);

xlim([-5,5])
ylim([-0.003,1.05])
title('Fig. 16.16: weight functions for a/symmetric outliers')
xlabel('$y$');ylabel('$w$');


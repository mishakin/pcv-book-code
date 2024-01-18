%% generate_step_corner
%
% [x,y] = generate_step_corner(N,sigma_e,sigma_n);
% 
% N         = number of points
% sigma_e   = process noise
% sigma_n   = observation noise
% Pout      = probability for outliers
% Mout      = maximal outlier
% type_out  = type otliers (0=symm., 1=asymm.)
% dens      = density of observed points in (0,1]
%

function [x,y,out_in,select,xs,ys] = ...
    generate_step_corner(Nd,sigma_e,sigma_n,Pout,Mout,type_out,dens,factor);

N=3*Nd;
x=zeros(3*Nd,1);
y=zeros(3*Nd,1);
h=ceil(Nd);

x1=h*ones(Nd,1);
x2=(1:Nd)';
x3=h-x2;
x=3+[x1;x2;x3]*factor;
%x=h*ones(N,1);
y=x+randn(N,1)*sigma_n;

out_in=zeros(N,1);
m=0;
for i=1:N
    if rand(1) < Pout
        out_in(i)=1;
        if type_out == 0
            y(i) = y(i) + 2*(rand(1)-0.5)*Mout;
        else
            y(i) = y(i) + rand(1)*Mout;
        end
    end
    
    if rand(1) < dens
        m=m+1;
        select(m)=i;
        xs(m)=x(i);
        ys(m)=y(i);
    end
end


return



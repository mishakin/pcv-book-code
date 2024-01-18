%% interpolate_bilinear


function z=interpolate_bilinear(S,x,y,dx,BB,sigma_h);

x=(x-BB(1))/dx+1;
y=(y-BB(2))/dx+1;

xr=floor(x);
yr=floor(y);
s=x-xr;
t=y-yr;
z=(1-s)*(1-t)*S(xr,yr)...
  +s*(1-t)*S(xr+1,yr)...
  +(1-s)*t*S(xr,yr+1)...
  +s*t*S(xr+1,yr+1); 
z=z*sigma_h;

[x-1,y-1,z,xr,yr,s,t,sqrt(S(xr,yr))*sigma_h];
return
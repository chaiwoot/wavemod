% 1D Acoustic modeling based on Second-Order wave Equation
% using the (2-2) Finite-Difference scheme with No Absorbing boundary
%
% dispersion condition: dx < c/f/npw, npw = number of points per wavelength
% stability condition: dt < dx/c
%
% Chaiwoot Boonyasiriwat, Mahidol University

clear;clc
npw = 20;
f = 100;
c = 1000;
dx_required = c/f/npw;
dx = dx_required;
dt_required = dx/c;
dt = dt_required;
cdt2 = c*c*dt*dt;
xmax = 100;
tmax = 1;
nx = floor(xmax/dx)+1;
nt = floor(tmax/dt)+1;
x = (0:nx-1)*dx;
t = (0:nt-1)*dt;
s = sin(2*pi*f*t);
p0 = zeros(nx,1);
p1 = p0;
p2 = p0;
alpha = (c*dt/dx)^2;
isx = round(nx/2);
for it=1:nt
    for ix=2:nx-1
        p2(ix) = (2-2*alpha)*p1(ix)-p0(ix)+alpha*(p1(ix+1)+p1(ix-1));
    end
    p2(isx) = p2(isx)+s(it)*cdt2;
    p0 = p1;
    p1 = p2;
    plot(x,p2);
    axis([0 xmax -1 1]);
    drawnow
end

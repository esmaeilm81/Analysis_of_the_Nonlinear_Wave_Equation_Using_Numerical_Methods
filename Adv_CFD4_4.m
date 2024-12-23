clc;
clear;
close all;

N=1001;
h=10/(N-1);
u_max=1;
xi=-4.995:0.01:4.995;
x=linspace(-5,5,1000);
a=length(xi);
u0=zeros(1,a);
u0(x>-1 & x<0)=-1;
u0(x>0 & x<1)=1;
t=0.5*h/u_max;
nt=floor(3/t);

u_p=u0;
f=@(u) 0.5*u^2
for n=1:nt
    u=u_p;
    for i=2:a-1
        flux_left=R_flux(u_p(i-1), u_p(i), f);
        flux_right=R_flux(u_p(i), u_p(i+1), f);
        u(i)=u_p(i)-(t/h)*(flux_right - flux_left);
    end
    u_p=u;
end
% Plot the result
plot(x, u, 'LineWidth', 2);
title('Solution at t = 2 s using Roe''s Method');
xlabel('x');
ylabel('u(x, t)');


function F=R_flux(uL,uR,f)
    F=(uL>uR)*(max(f(uL),f(uR)))+...
    (uL<=uR)*(min(f(uL),f(uR)));
end
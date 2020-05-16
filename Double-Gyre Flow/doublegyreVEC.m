

function [dy] = doublegyreVEC(t,yin,A,eps,om)
%yin is a 3D matrix: thr first coloum is 2 (veclocity components)
% the second and third coloums are the mesh dimensions. e.g., yin = (2,100,50)
x = yin(1,:,:);
y = yin(2,:,:);
%define streamfunction parameters (from Shadden 2005 Physica D)
% A =  0.1; 
%eps = 0.25; 
%om = 2*pi/10;  
%define the actual velovities
u = zeros(size(x)); 
v = u;
%define frequency parameters
a = eps * sin(om * t);
b = 1 - 2 * a;
%defne the time varying frequency
f = a * x.^2 + b * x;
df = 2 * a * x + b;
%compute the velocities
u = -pi * A * sin(pi * f) .* cos(pi * y);
v =  pi * A * cos(pi * f) .* sin(pi * y) .* df;

dy = [u;v];
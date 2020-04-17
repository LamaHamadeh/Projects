close all;
clear all;


%define x range
x = -1:0.01:1; 

%defne functions on x domain
u = exp(x).*sin(5*x); %1st function
ux = exp(x).*sin(5*x)+5*exp(x).*cos(5*x); %1st derivative
uxx = -24*exp(x).*sin(5*x)+10*exp(x).*cos(5*x); %2nd derivative

v = sech(x); %2nd function
vx = -sech(x).*tanh(x); %1st derivative
vxx = sech(x)-2*sech(x).^3; %2nd derivative

%------                                    ------%
% How good chebychev matrix do to real derivatives? 
%------                                    ------%

%define the number of points
N=40;
%apply the Chebychev differentiation matrix
[D,x2] = Cheb_Diff_Matrix(N); %x2 is the new domain that have N points

%define functions on the points N
u2 = exp(x2).*sin(5*x2); %1st function
u2x = D*u2; %1st derivative
u2xx = (D^2)*u2; %2nd derivative

v2 = sech(x2); %2nd function
v2x = D*v2; %1st derivative
v2xx = (D^2)*v2; %2nd derivative

%plotting
%----------
figure(1) %for the first function
plot(x,u,'r',x2,u2,'ro',x,ux,'m',x2,u2x,'mo',x,uxx,'k',x2,u2xx,'ko')
figure(2) %for the second function
plot(x,v,'r',x2,v2,'ro',x,vx,'m',x2,v2x,'mo',x,vxx,'k',x2,v2xx,'ko')


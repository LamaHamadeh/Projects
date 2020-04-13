
clear all; 
close all; 

                    %%%Analytical Solution

%Defining temperatures
%body
Tbody_C = 37; %celsius
Tbody_K = Tbody_C + 273.15; %kelvin
%coffee temp1
Tcoffee1_C = 60; %celsius
Tcoffee1_K = Tcoffee1_C + 273.15; %kelvin
% %coffee temp2
% Tcoffee2_C = 82; %celsius
% Tcoffee2_K = Tcoffee2_C + 273.15; %kelvin

% %space variable
% L=2; %mm
% dx=0.05;
% x = 0:dx:L;
% 
% %time variable
% tmax=4; %sec
% dt=0.01;
% t = 0:dt:tmax;
% 
% %themral conductivity
% k1 = 0.209; %W/mK for epidermis
% %k2 = 0.322; %W/mK for dermis
% 
% %analytical solution
% [X,T] = meshgrid(x,t);
% U0 = ((Tcoffee1_K-Tbody_K)*X)/L+Tbody_K;
% tic
% syms n
% U = U0 + symsum(2*(Tbody_K-Tcoffee1_K)*((-1)^n)/pi*n.* sin(n*pi*X/L)...
%     .*exp(-k1*T*(n*pi/L)^2),n,1,1000);
% toc
%U = zeros(size(X,1),size(T,2));

% for n = 1:1000
%      Bn = 2*(Tbody_K-Tcoffee1_K)*((-1)^n)/pi*n;
%      %Bn = -2*(-1)^n/(n*pi)^2;
%      Un = Bn .* sin(n*pi*X/L).*exp(-k1*T*(n*pi/L)^2);
%      U = U0 + Un;
%      U = U - 273.15;
% end

%plotting
%convert to celsius

%figure(2)
%plot(X,U(1,:),'b',X,U(100,:),'r',X,U(200,:),'g',X,U(300,:),'k',X,U(400,:),'m')


                     %----------------
                    %%%Numerical Solution

%initialise the grid
%space variable
L = 2;
n=400;
x2 = linspace(0,L,n+1);
x = x2(1:n);
dx = 0.5;
%time variable
dt = 0.02;
Time = 4;
t = 0:dt:Time;
%CFL to ensure stability
CFL = dt/dx^2;

%define solution matrix
Tsol = zeros(length(x),length(t));
%initial condition
Tsol(:,1)=Tbody_K;
T0 = Tsol(:,1);
%boundary conditions
Tsol(1,:) = Tbody_K;
Tsol(n,:) = Tcoffee1_K;

%laplacian matrix
e1 = ones(n,1);
A = spdiags([e1 -2*e1 e1],[-1 0 1],n,n);
A(1,n) = 1;
A(n,1) = 1;

%solve using Forward Euler method
for j = 1:length(t)-1
    T1 = T0+CFL*A*T0;
    T0 = T1;
    Tsol(:,j+1) = T0; 
end

%back to celsius
Tsol = Tsol-273.15;

plot(x,Tsol(:,200))
%plotting
waterfall(x,t,Tsol');
map = [0 0 0]; %BW
colormap(map);
 







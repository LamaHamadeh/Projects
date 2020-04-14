
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


                    %%%Numerical Solution
%initialise numerical variables
%space 
L = 2; %skin depth
dx = 0.2; %space step
nx = L/dx; %number of space points
x = linspace(0,L,nx); %space variable
%time 
Time = 4; %maximum time 
dt = 0.02; %time step
nt = Time/dt; %number of time points
t = linspace(0,Time,nt); %time variable

%themral conductivity
k = 0.209; %epidermis dominated effect
%k = 0.322; %dermis dominated effect

%CFL number to ensure stability (less or equal to 0.5)
CFL = dt/dx^2;

%Define boundary temperatures
%body temperature 
Tbody_C = 37; %in celsius
Tbody_K = Tbody_C + 273.15; %kelvin
%coffee temperature 1
% Tcoffee_C = 60; %in celsius
% Tcoffee_K = Tcoffee_C + 273.15; %kelvin
%coffee temperature 2
Tcoffee_C = 82; %in celsius
Tcoffee_K = Tcoffee_C + 273.15; %kelvin

%initial condition
Tini = zeros(nx,1);
Tini(:,1) = Tbody_K;

%boundary conditions
%BC at x=0
T0 = zeros(1,nt); %define a vector
T0(1,:) = Tbody_K; 
%BC at x=L
TL = zeros(1,nt); %define a vector
TL(1,:) = Tcoffee_K;

%construct the solution matrix
Tsol = zeros(nx,nt); %define the dimensions of the solution matrix
Tsol(1,:) = T0; %insert the first BC
Tsol(nx,:) = TL; %insert the second BC
Tsol(:,1) = Tini; %insert the IC

%set up the iteration process for forward Euler method
for i = 1:nt-1
    for j = 2:nx-1
        Tsol(j,i+1) = Tsol(j,i) + k * CFL *(Tsol(j+1,i)-2*Tsol(j,i)+Tsol(j-1,i));
    end
end

%convert the solution to celsius
Tsol = Tsol-273.15;

%plotting
%for different time steps
figure(1)
plot(x,Tsol(:,1),'b',x,Tsol(:,50),'r',x,Tsol(:,100),'g',...
    x,Tsol(:,150),'m',x,Tsol(:,200),'k','LineWidth',2)
xlabel('Depth (mm)')
ylabel('Temperature (C^0)')
title(['Thermal conductivity= ', num2str(k), ' and coffee temperature= ', ...
    num2str(Tcoffee_C),' C^0'])
legend('0 sec','1 sec','2 sec','3 sec','4 sec','Location','NorthWest')
set(gca,'FontSize',16)
axis([0 2 37 85]);

%3D plot
figure(2)    
mesh(x,t,Tsol')
colormap(jet)
shading interp
xlabel('Depth (mm)')
ylabel('Time (Sec)')
zlabel('Temperature (C^0)')
colorbar
set(gca,'FontSize',16)








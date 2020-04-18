
clear all;
close all;

                    %%%Analytical Solution

%initialise spacetime variables
%space 
L = 2; %skin depth
dx = 0.2   ; %space step
nx = L/dx; %number of space points
x = 0:dx:L; %space variable
%time 
Time = 4; %maximum time 
dt = 0.02; %time step
nt = Time/dt; %number of time points
t = 0:dt:Time; %time variable

%themral conductivity
%k = 0.209; %epidermis dominated effect
k = 0.322; %dermis dominated effect

%Define boundary temperatures
%body temperature 
Tbody_C = 37; %in celsius 
Tbody_K = Tbody_C + 273.15; %kelvin
%coffee temperature 1
Tcoffee_C = 60; %in celsius
Tcoffee_K = Tcoffee_C + 273.15; %kelvin
%coffee temperature 2
% Tcoffee_C = 82; %in celsius
% Tcoffee_K = Tcoffee_C + 273.15; %kelvin

%Define anonymous functions
C_n = @(Tbody_K,Tcoffee_K,n) 2*(Tcoffee_K-Tbody_K)*((-1)^n)/(pi*n);
C_0 = @(Tbody_K,Tcoffee_K,x) (Tcoffee_K-Tbody_K)*x/L+Tbody_K;
Tn = @(n,x,t,k,Tbody_K,Tcoffee_K) C_n(Tbody_K,Tcoffee_K,n).*sin(n*pi*x/L).*exp(-k*t*(n*pi/(L))^2);

%Define a mesh for space and time
[X,T] = meshgrid(x,t);

%initialise the sum to zero
Tsum = zeros(size(C_0(Tcoffee_K,Tbody_K,X)));
%use definite loop instead of indefinite sum
for i=1:1000
    Tsum = Tsum + Tn(i,X,T,k,Tbody_K,Tcoffee_K);%computing the sum
end

%compute the final solution
Tsol = C_0(Tbody_K,Tcoffee_K,X)+Tsum;

%convert the solution to celsius
Tsol = (Tsol-273.15)';

%plotting
%for different time steps
figure(1)
plot(x,Tsol(:,1),'b',x,Tsol(:,50),'r',x,Tsol(:,150),'g',...
    x,Tsol(:,175),'m',x,Tsol(:,200),'k','LineWidth',2)
xlabel('Depth (mm)')
ylabel('Temperature (C^0)')
title(['Thermal conductivity= ', num2str(k), ' and coffee temperature= ', ...
    num2str(Tcoffee_C),' C^0'])
legend('0 sec','1 sec','2 sec','3 sec','4 sec','Location','NorthWest')
set(gca,'FontSize',16)

figure(2)
plot(t,Tsol(1,:),'b',t,Tsol(3,:),'r',t,Tsol(6,:),'g',...
    t,Tsol(8,:),'m',t,Tsol(10,:),'k','LineWidth',2)
xlabel('Time (sec)')
ylabel('Temperature (C^0)')
title(['Thermal conductivity= ', num2str(k), ' and coffee temperature= ', ...
    num2str(Tcoffee_C),' C^0'])
legend('0 mm','0.5 mm','1.5 mm','1.75 mm','2 mm','Location','NorthWest')
set(gca,'FontSize',16)


%                     %%%Numerical Solution
% %initialise numerical variables
% %space 
% L = 2; %skin depth
% dx = 0.2; %space step
% nx = L/dx; %number of space points
% x = linspace(0,L,nx); %space variable
% %time 
% Time = 4; %maximum time 
% dt = 0.02; %time step
% nt = Time/dt; %number of time points
% t = linspace(0,Time,nt); %time variable
% 
% %themral conductivity
% %k = 0.209; %epidermis dominated effect
% k = 0.322; %dermis dominated effect
% 
% %CFL number to ensure stability (less or equal to 0.5)
% CFL = dt/dx^2;
% 
% %Define boundary temperatures
% %body temperature 
% Tbody_C = 37; %in celsius 
% Tbody_K = Tbody_C + 273.15; %kelvin
% %coffee temperature 1
% Tcoffee_C = 60; %in celsius
% Tcoffee_K = Tcoffee_C + 273.15; %kelvin
% %coffee temperature 2
% % Tcoffee_C = 82; %in celsius
% % Tcoffee_K = Tcoffee_C + 273.15; %kelvin
% 
% %initial condition
% Tini = zeros(nx,1);
% Tini(:,1) = Tbody_K;
% 
% %boundary conditions
% %BC at x=0
% T0 = zeros(1,nt); %define a vector
% T0(1,:) = Tbody_K; 
% %BC at x=L
% TL = zeros(1,nt); %define a vector
% TL(1,:) = Tcoffee_K;
% 
% %construct the solution matrix
% Tsol = zeros(nx,nt); %define the dimensions of the solution matrix
% Tsol(1,:) = T0; %insert the first BC
% Tsol(nx,:) = TL; %insert the second BC
% Tsol(:,1) = Tini; %insert the IC
% 
% %set up the iteration process for forward Euler method
% for j = 1:nt-1 
%     for i = 2:nx-1
%         Tsol(i,j+1) = Tsol(i,j) + k * CFL *(Tsol(i+1,j)-2*Tsol(i,j)+Tsol(i-1,j));
%     end
% end
% 
% %convert the solution to celsius
% Tsol = Tsol-273.15;
% 
% %plotting
% %for different time steps
% figure(1)
% plot(x,Tsol(:,1),'b',x,Tsol(:,50),'r',x,Tsol(:,100),'g',...
%     x,Tsol(:,150),'m',x,Tsol(:,200),'k','LineWidth',2)
% xlabel('Depth (mm)')
% ylabel('Temperature (C^0)')
% title(['Thermal conductivity= ', num2str(k), ' and coffee temperature= ', ...
%     num2str(Tcoffee_C),' C^0'])
% legend('0 sec','1 sec','2 sec','3 sec','4 sec','Location','NorthWest')
% set(gca,'FontSize',16)
% 
% figure(2)
% plot(t,Tsol(1,:),'b',t,Tsol(3,:),'r',t,Tsol(6,:),'g',...
%     t,Tsol(8,:),'m',t,Tsol(10,:),'k','LineWidth',2)
% xlabel('Time (sec)')
% ylabel('Temperature (C^0)')
% title(['Thermal conductivity= ', num2str(k), ' and coffee temperature= ', ...
%     num2str(Tcoffee_C),' C^0'])
% legend('0 mm','0.5 mm','1.5 mm','1.75 mm','2 mm','Location','NorthWest')
% set(gca,'FontSize',16)
% 
% 



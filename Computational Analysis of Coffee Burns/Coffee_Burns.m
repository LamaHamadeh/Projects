
clear all; 
close all; 

                    %%%Analytical Solution



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
%k = 0.209; %epidermis dominated effect
k = 0.322; %dermis dominated effect

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

figure(2)
plot(t,Tsol(1,:),'b',t,Tsol(3,:),'r',t,Tsol(6,:),'g',...
    t,Tsol(8,:),'m',t,Tsol(10,:),'k','LineWidth',2)
xlabel('Time (sec)')
ylabel('Temperature (C^0)')
title(['Thermal conductivity= ', num2str(k), ' and coffee temperature= ', ...
    num2str(Tcoffee_C),' C^0'])
legend('0 mm','0.5 mm','1.5 mm','1.75 mm','2 mm','Location','NorthWest')
set(gca,'FontSize',16)







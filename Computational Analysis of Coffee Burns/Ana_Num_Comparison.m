
clear all;
close all;

                    %%%Define Global Parameters
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
k = 0.209; %epidermis dominated effect
%k = 0.322; %dermis dominated effect

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


                    %%%Analytical Solution

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
Tsol_ana = C_0(Tbody_K,Tcoffee_K,X)+Tsum;

%convert the solution to celsius
Tsol_ana = (Tsol_ana-273.15)';

                    %%%Numerical Solution

%CFL number to ensure stability (less or equal to 0.5)
CFL = dt/dx^2;

%initial condition
Tini = zeros(nx+1,1);
Tini(:,1) = Tbody_K;

%boundary conditions
%BC at x=0
T0 = zeros(1,nt+1); %define a vector
T0(1,:) = Tbody_K; 
%BC at x=L
TL = zeros(1,nt+1); %define a vector
TL(1,:) = Tcoffee_K;

%construct the solution matrix
Tsol_num = zeros(nx+1,nt+1); %define the dimensions of the solution matrix
Tsol_num(1,:) = T0; %insert the first BC
Tsol_num(nx+1,:) = TL; %insert the second BC
Tsol_num(:,1) = Tini; %insert the IC

%set up the iteration process for forward Euler method
for j = 1:nt
    for i = 2:nx
        Tsol_num(i,j+1) = Tsol_num(i,j) + k * CFL *(Tsol_num(i+1,j)-2*Tsol_num(i,j)+Tsol_num(i-1,j));
    end
end

%convert the solution to celsius
Tsol_num = Tsol_num-273.15;


                        %%% Global (cumulative) Error Estimation of the Numerical Scheme 
                        %%% Global Disretisation Error (Error = Ana(i) - num(i)) where i is the time step

%Note: if we denote:

%A to be the analytical solution of the partial differential equations
%D the exact solution of the difference equation
%N the numerical solution from a computer, i.e., iterative solution

%Then we can define:

%Discretisation error = A - D
%Round-off error = N - D
%Global error of a numerical scheme = A - N

%Root mean square error
for i = 1:length(t)
    Error = sqrt(sum((Tsol_ana(:,i) - Tsol_num(:,i)).^2)/numel(x));
end

                        %%%Plotting
                        
%for different time steps
figure(1)
plot(x,Tsol_ana(:,1),'b',x,Tsol_num(:,1),'bo',...
    x,Tsol_ana(:,50),'r',x,Tsol_num(:,50),'ro',...
    x,Tsol_ana(:,100),'g',x,Tsol_num(:,100),'go',...
    x,Tsol_ana(:,150),'m',x,Tsol_num(:,150),'mo',...
    x,Tsol_ana(:,200),'k',x,Tsol_num(:,200),'ko','LineWidth',2)
xlabel('Depth (mm)')
ylabel('Temperature (C^0)')
title(['Thermal conductivity= ', num2str(k), ' and coffee temperature= ', ...
    num2str(Tcoffee_C),' C^0'])
legend('0 sec Analytical','0 sec Numerical',...
    '1 sec Analytical','1 sec Numerical',...
    '2 sec Analytical','2 sec Numerical',...
    '3 sec Analytical','3 sec Numerical',...
    '4 sec Analytical','4 sec Numerical',...
    'Location','NorthWest')
set(gca,'FontSize',16)
%for different spatial depths
figure(2)
plot(t,Tsol_ana(1,:),'b',t,Tsol_num(1,:),'bo',...
    t,Tsol_ana(6,:),'r',t,Tsol_num(6,:),'ro',...
    t,Tsol_ana(10,:),'g',t,Tsol_num(10,:),'go','LineWidth',2)
xlabel('Time (sec)')
ylabel('Temperature (C^0)')
title(['Thermal conductivity= ', num2str(k), ' and coffee temperature= ', ...
    num2str(Tcoffee_C),' C^0'])
legend('0 mm Analytical','0 mm Numerical',...
    '0.5 mm Analytical','0.5 mm Numerical',...
    '2 mm Analytical','2 mm Numerical',...
    'Location','NorthWest')
set(gca,'FontSize',16)
% Plot imagesc plots
figure(3)
%analytical solution
subplot(1,2,1)
imagesc(Tsol_ana)
colorbar;
axis square
title('Analytical')
xlabel('Time')
ylabel('Depth')
set(gca,'FontSize',14)
colormap summer
%numerical solution
subplot(1,2,2)
imagesc(Tsol_num)
colorbar;
axis square 
title('Numerical')
xlabel('Time')
ylabel('Depth')
set(gca,'FontSize',14)
colormap summer


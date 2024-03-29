clear all; 
close all;

%One-way wave equation
%---------------------

%initialise the grid, time,space and CFL number
%time
Time = 2;
dt = 0.005; %makes CFL =2 and both schemes to be instable
%dt = 0.1; %makes CFL = 1 and Leap Frog to be stable. Forward Euler Scheme
%is unstabhle for any value of CFL
t = 0:dt:Time;

%space
L = 20;
n = 150;
x2 = linspace(-L/2,L/2,n+1);
x = x2(1:n);
dx = x(2)-x(1);

%CFL number (key parameter)
CFL_wave = dt/dx;

%NOTE: normally, the smaller dx and dt are, the m ore your system is converged.
%to test your convergence, always try to chop yout dt in half (or increase the max time). If you get 
%the same solution (same sybnamics), then your system is correctly converged.
%Otherwise, your system is instable.

                                %--------------------

%Forward Euler Scheme
%---------------------
%initial condition
u0 = exp(-x.^2).'; %the transpose is to set the values as a column
usol(:,1) = u0; %to save the result in the first column
%plot(x,u0);

%u_x matrix (derivative matrix)
e1 = ones(n,1);
A = spdiags([-e1 e1], [-1 1], n, n);
A(1,n) = -1;
A(n,1) = 1;

%set up the iteration process 
for j = 1:length(t)-1
    u1 = u0 + (CFL_wave/2)*A*u0; %future time point
    u0 = u1;   
    usol(:,j+1) = u1;
end

%plot the one way wave
figure(1)
waterfall(x,t,usol');
map = [ 0 0 0];
colormap(map);
xlabel('Distance')
ylabel('Time')
zlabel('Wave Solution')
% surf(x,t,usol');
% colormap(winter)

%this scheme is unstable whatever dt or dx values are. this should be clear
%in the plotted solution.
                                %--------------------

%Leap frog Scheme (2,2)
%----------------------
%initial condition
%this is a two-step procedure therefore we nedd two initial conditions.
u0 = exp(-x.^2).'; %the transpose is to set the values as a column
u1 = exp(-(x+dt).^2).';
usol(:,1) = u0; 
usol(:,2) = u1;
%plot(x,u0);
 
%u_x matrix
e1 = ones(n,1);
A = spdiags([-e1 e1], [-1 1], n, n);
A(1,n) = -1;
A(n,1) = 1;
 
%set up the iteration process 
for j = 1:length(t)-2
    u2 = u0 + CFL_wave*A*u1;
    u0 = u1; %my past is my current
    u1 = u2; %my current is my future
    usol(:,j+2) = u2;
end
% 
%plot the one way wave
figure(2)
waterfall(x,t,usol');
map = [ 0 0 0];
colormap(map);
xlabel('Distance')
ylabel('Time')
zlabel('Wave Solution')

%CFL here in this scheme should be CFL<=1.
                                %--------------------

%Heat equation
%--------------

CFL_heat = dt/dx^2;

%Forward Euler Scheme
%---------------------
%initial condition
u0 = exp(-x.^2).'; %the transpose is to set the values as a column
u1 = exp(-(x).^2).';
usol(:,1) = u0; 
usol(:,2) = u1;

%u_x matrix (derivative matrix)
e1 = ones(n,1);
A = spdiags([e1 -2*e1 e1], [-1 0 1], n, n);
A(1,n) = 1;
A(n,1) = 1;

%set up the iteration process 
for j = 1:length(t)-1
    u1 = u0 + CFL_heat*A*u0; %future time point
    u0 = u1;   
    usol(:,j+1) = u1;
end

%plot the one way wave
figure(3)
waterfall(x,t,usol');
map = [ 0 0 0];
colormap(map);
xlabel('Distance')
ylabel('Time')
zlabel('Wave Solution')
% surf(x,t,usol');
% colormap(winter)


                                %--------------------


%Leap Frog Scheme (2,2)
%----------------------

%initial condition
%this is a two-step procedure therefore we nedd two initial conditions.
u0 = exp(-x.^2).'; %the transpose is to set the values as a column
u1 = exp(-(x).^2).';
usol(:,1) = u0; 
usol(:,2) = u1;
%plot(x,u0);
 
%u_x matrix
e1 = ones(n,1);
A = spdiags([e1 -2*e1 e1], [-1 0 1], n, n);
A(1,n) = 1;
A(n,1) = 1;
 
%set up the iteration process 
for j = 1:length(t)-2
    u2 = u0 +2* CFL_heat*A*u1;
    u0 = u1; %my past is my current
    u1 = u2; %my current is my future
    usol(:,j+2) = u2;
end
% 
%plot the one way wave
figure(4)
waterfall(x,t,usol');
map = [ 0 0 0];
colormap(map);
xlabel('Distance')
ylabel('Time')
zlabel('Wave Solution')

%Leap frog in the second derivative is always not stable

                                %--------------------


%Non-linear Schodinger Equation
%-------------------------------

%Leap Frog Scheme (2,2)
%----------------------

%initial condition
%this is a two-step procedure therefore we nedd two initial conditions.
u0 = exp(-x.^2).'; %the transpose is to set the values as a column
u1 = exp(-(x).^2).';
usol(:,1) = u0; 
usol(:,2) = u1;
%plot(x,u0);
 
%u_x matrix
e1 = ones(n,1);
A = spdiags([e1 -2*e1 e1], [-1 0 1], n, n);
A(1,n) = 1;
A(n,1) = 1;
 
%set up the iteration process 
for j = 1:length(t)-2
    u2 = u0 -(1i*CFL_heat*A*u1-1i*2*dt*(conj(u1).*u1).*u1);
    u0 = u1; %my past is my current
    u1 = u2; %my current is my future
    usol(:,j+2) = u2;
end
% 
%plot the one way wave
figure(5)
waterfall(x,t,abs(usol'));
map = [ 0 0 0];
colormap(map);
xlabel('Distance')
ylabel('Time')
zlabel('Wave Solution')
















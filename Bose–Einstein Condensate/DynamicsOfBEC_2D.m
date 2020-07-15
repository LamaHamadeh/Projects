%%Still Working on it!


clear all;
close all;

%Spatial variable on x direction
Lx=5; %domain length
delta=0.05; %spatial step size
xmin=-Lx; %minimum boundary
xmax=Lx; %maximum boundary 
Nx=(xmax-xmin)/delta; %number of spatial points
x=linspace(xmin,xmax,Nx); %spatial vector
%--------------

%Spatial variable on y direction
Ly=5; %domain length
delta=0.05; %spatial step size
ymin=-Ly; %minimum boundary
ymax=Ly; %maximum boundary M
Ny=(ymax-ymin)/delta; %number of spatial points
y=linspace(ymin,ymax,Ny); %spatial vector
%--------------

%Total matrix size
N = Nx * Ny;
%--------------

%Create a meshgrid to let Matlab know the x and y directions
[X,Y] = meshgrid(x,y);
%--------------

%Define sinusoidal trapping potential
%parameters
A1 = 0.5;
B1 = -0.5;
A2 = 0.5;
B2 = -0.5;
%potential
V = (A1.*(sin(X)).^2+B1).*(A2.*(sin(Y)).^2+B2)...
    +(sin(X)).^2+(sin(Y)).^2;
%plotting trapping potential
% get the corners of the domain in which the data occurs.
min_x = min(min(X));
min_y = min(min(Y));
max_x = max(max(X));
max_y = max(max(Y));
% the image data you want to show as a plane.
planeimg = V;
% set hold on so we can show multiple plots / surfs in the figure.
figure; 
hold on;
% do a normal surface plot.
surf(X,Y,V);
shading interp
% desired z position of the image plane.
imgzposition = -2;
% plot the image plane using surf.
surf([min_x max_x],[min_y max_y],repmat(imgzposition, [2 2]),...
    planeimg,'facecolor','texture')
% set a colormap for the figure.
colormap default;
% set the view angle.
view(45,30);
%labelling
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$V(x,y)$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
%zlim([-3 2])
%--------------

%Time variable
dt=0.01;  %time step size
tmin=0; %minimum time
tmax=1; %maximum time
nt=(tmax-tmin)/dt; %number of time points
tspan=linspace(tmin,tmax,nt); %time variable
%--------------



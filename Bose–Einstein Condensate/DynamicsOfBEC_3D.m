
clear all;
close all;

%Spatial variable on x direction
Lx=3; %domain length
delta=0.05; %spatial step size
xmin=-Lx; %minimum boundary
xmax=Lx; %maximum boundary 
Nx=(xmax-xmin)/delta; %number of spatial points
x=linspace(xmin,xmax,Nx); %spatial vector
%--------------

%Spatial variable on y direction
Ly=3; %domain length
delta=0.05; %spatial step size
ymin=-Ly; %minimum boundary
ymax=Ly; %maximum boundary M
Ny=(ymax-ymin)/delta; %number of spatial points
y=linspace(ymin,ymax,Ny); %spatial vector
%--------------

%Spatial variable on z direction
Lz=3; %domain length
delta=0.05; %spatial step size
zmin=-Lz; %minimum boundary
zmax=Lz; %maximum boundary M
Nz=(zmax-zmin)/delta; %number of spatial points
z=linspace(zmin,zmax,Nz); %spatial vector
%--------------

%Total matrix size
N = Nx * Ny * Nz;
%--------------

%Create a meshgrid to let Matlab know the x and y directions
[X,Y,Z] = meshgrid(x,y,z);
%--------------

%Define sinusoidal trapping potential
%parameters
A1 = 0.5;
B1 = -0.5;
A2 = 0.5;
B2 = -0.5;
A3 = 0.5;
B3 = -0.5;
%potential
V = (A1.*(sin(X)).^2+B1).*(A2.*(sin(Y)).^2+B2).*(A3.*(sin(Z)).^2+B3)...
    +(sin(X)).^2+(sin(Y)).^2+(sin(Z)).^2;
figure;
%using slice
% min_x = min(min(X));
% min_y = min(min(Y));
% min_z = min(min(Z));
% max_x = max(max(X));
% max_y = max(max(Y));
% max_z = max(max(Z));
% xslice = [min_x max_x];   
% yslice = [min_y max_y];
% zslice = [min_z max_z];
% slice(X,Y,Z,V,xslice,yslice,zslice,'nearest')
% shading interp
% xlabel('$x$','Interpreter','latex')
% ylabel('$y$','Interpreter','latex')
% zlabel('$z$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'FontSize',16)

%using isosurface
isovalue = 2.9;
isosurface(X,Y,Z,V,isovalue)
view(3); 
axis tight
camlight 
lighting gouraud
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
%--------------

% %Time variable
% dt=0.01;  %time step size
% tmin=0; %minimum time
% tmax=1; %maximum time
% nt=(tmax-tmin)/dt; %number of time points
% tspan=linspace(tmin,tmax,nt); %time variable
% %--------------



%Date: 28/01/2020
%----------------

%Description: 
%------------
%This file solve the partial differential equation of the time-dependent 
%Thermal distribution associated with one-dimensional heat transfer in a 
%stationary, homogeneous medium using Finite Difference Method (FDM).


clear; close all;

%Defining Numerical parameters:
%--------------------------------
%Spatial variable on x direction
%--------------------------------
Lx=1; %LED's length
delta=0.05; %spatial step size
xmin=-Lx; %minimum boundary
xmax=Lx; %maximum boundary 
Nx=(xmax-xmin)/delta; %number of spatial points
x=linspace(xmin,xmax,Nx); %spatial vector

%Spatial variable on y direction
%--------------------------------
Ly=1; %LED's length
delta=0.05; %spatial step size
ymin=-Ly; %minimum boundary
ymax=Ly; %maximum boundary M
Ny=(ymax-ymin)/delta; %number of spatial points
y=linspace(ymin,ymax,Ny); %spatial vector

%Total matrix size
%-----------------
N = (Nx * Ny);

%Time variable
%-------------
dt=0.01;  %time step size
tmin=0; %minimum time
tmax=1; %maximum time
nt=(tmax-tmin)/dt; %number of time points
tspan=linspace(tmin,tmax,nt); %time variable

%Create a meshgrid to let Matlab know the x and y directions
%-----------------------------------------------------------
[X,Y] = meshgrid(x,y);

% Defining initial state:
%-------------------------
%2D Gaussian Function 
T0=exp(-(X.^2+Y.^2)); 
%figure;
%imagesc(T0)

%2D square function
% T0 = (2/pi)*[(atan(exp(-1i*(pi/2).*X))+atan(exp(1i*(pi/2).*X)))+ ...
%     (atan(exp(-1i*(pi/2).*Y))+atan(exp(1i*(pi/2).*Y)))];

%reshape the initial condition to a vector
%------------------------------------------
T_reshape = reshape(T0,N,1); 

% Constructing the spatial matrix using Finite-Difference Method (FDM)
%------------------------------------------------------------------------
A=zeros(Nx,Nx); %matrix of zero elements
I = eye(Nx);

%the diagonal elements
for m=1:Nx %the number of rows
    for n=1:Nx  %the number of colomns
if (m==n)
    A(m,n)=-2/delta^2;  %the value of each diagonal element 
end
%Boundary conditions: A(1,N)==A(N,1)
if(n==Nx)&&(m==1)
    A(m,n)=1;
end

if(n==1)&&(m==Nx)
    A(m,n)=1;
end
    end
end
%the off-diagonal elements
for n=1:Nx-1       
    A(n+1,n)=1/delta^2; %the value of each lower off-diagonal elements
end 
for n=2:Nx        
    A(n-1,n)=1/delta^2; %the value of each upper off-diagonal element
end

%sparse A: extracts the nonzero diagonals from m-by-n matrix A and returns 
%them as the columns in min(m,n)-by-p matrix B,
% A = spdiags(A);

%create the 2D matrix
%---------------------
B = kron(A,I)+kron(I,A);

% Solve the one dimensional, time-dependent and partial differential equation
%---------------------------------------------------------------------------
[Time,Tem]=ode45('dTDistribution',tspan,T_reshape,[],B,delta); 
%[] refers to tolerence
%the rows of Tem are the solutions for each given time


%reshape the solution from a vector to a 2D x-y plane and plot it
%-----------------------------------------------------------------
h = figure;
axis tight manual
filename = 'FDM2D_gaussian.gif';

for k = 1:5
    for i =1:size(Tem,1)
        %reshape the solution and plot it at each time step
        TFinal = reshape(Tem(i,:),Nx,Ny);
        
        %draw surf plots
        surf(X,Y,TFinal), shading interp, colormap hot, colorbar
        %axis labeling
        xlabel('x')
        ylabel('y')
        zlabel('Diffusion')
        %axis limits
        xlim([-1.5 1.5]) 
        ylim([-1.5 1.5])
        zlim ([0, 1]); %This is very important so we can lock the zoom-in behaviour
        %grid
        grid off
        
        %capture the plot as an image
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        %write to the GIF file
        if k == 1
            imwrite(imind,cm,filename, 'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename, 'gif','WriteMode','append');
        end
    end
end

%Stability analysis
%-------------------
c = 1; %diffusion coefficient
CFL = c*dt/delta;
%check stability
if CFL >=1
    fprintf('The solution is not stable, Try again with different time and space discretisation steps.')
end
if CFL < 1
    fprintf('Stable solution!')
end















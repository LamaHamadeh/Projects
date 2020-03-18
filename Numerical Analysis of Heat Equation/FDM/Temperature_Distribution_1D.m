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
%Spatial variable
%-----------------
L=1; %LED's length
dx=0.01; %spatial step size
xmin=-L; %minimum boundary
xmax=L; %maximum boundary 
N=(xmax-xmin)/dx; %number of spatial points
x=linspace(xmin,xmax,N); %spatial vector


%Time variable
%-------------
dt=0.005;  %time step size
tmin=0; %minimum time
tmax=1; %maximum time
nt=(tmax-tmin)/dt; %number of time points
tspan=linspace(tmin,tmax,nt); %time variable


% Defining initial state:
%-------------------------
%Gaussian Function 
%T0=exp(-(x.^2));
%Square Function
T0 = (2/pi)*(atan(exp(-1i*(pi/2).*x))+atan(exp(1i*(pi/2).*x))); 

% Constructing the spatial matrix using Finite-Difference Method (FDM)
%------------------------------------------------------------------------
A=zeros(N,N); %matrix of zero elements
%the diagonal elements
for m=1:N %the number of rows
    for n=1:N  %the number of colomns
if (m==n)
    A(m,n)=-2/dx^2;  %the value of each diagonal element 
end
%Boundary conditions: A(1,N)==A(N,1)
if(n==N)&&(m==1)
    A(m,n)=1;
end

if(n==1)&&(m==N)
    A(m,n)=1;
end
    end
end
%the off-diagonal elements
for n=1:N-1       
    A(n+1,n)=1/dx^2; %the value of each lower off-diagonal element
end 
for n=2:N        
    A(n-1,n)=1/dx^2; %the value of each upper off-diagonal element
end

%if A is a big matrix, we can sparse it
% A = spdiags(A);

% Solve the one dimensional, time-dependent and partial differential equation
%---------------------------------------------------------------------------
[Time,Tem]=ode45('dTDistribution',tspan,T0,[],A,dx); 
%[] refers to tolerence
%the rows of Tem are the solutions for each given time

% Plotting solutions
%---------------------
%2D plots
for i = 1:length(Tem)
    figure(1)
    plot(x,Tem(i,:))
    hold on
end

%contour plot
figure(2)
contour(x,tspan,Tem,200)



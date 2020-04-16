

clear all;
close all;

           %%% Define Discretised spatial and time variables %%% 

%Space
%-----
%Spatial variable on x direction
Nx=64; %number of spatial points
x2 = linspace(-10,10,Nx+1); %x variable
x = x2(1:Nx);
Lx = 20;
%Spatial variable on y direction
Ny=64; %number of spatial points
y2 = linspace(-10,10,Ny+1); %y variable
y = y2(1:Ny);
Ly = 20;
delta=20/Nx; %spatial step size
%Total matrix size
N = Nx*Ny;
%Time 
%-----
dt=3;  %time step size
tspan=0:dt:47; %time variable
%Viscosity
nu = 0.001;

                    %%% Wave vector disretisation %%%
                   
%x direction
kx = (2*pi/Lx)*[0:(Nx/2-1) (-Nx/2):-1]'; 
kx(1) = 10^(-6);
%y direction
ky = (2*pi/Ly)*[0:(Ny/2-1) (-Ny/2):-1]'; 
ky(1) = 10^(-6);
%to give kx and ky the sense of direction
[KX,KY] = meshgrid(kx,ky);
Kderv = KX.^2+KY.^2;

                	     %%% Initial state %%%

%initial consitions
[X,Y] = meshgrid(x,y);
%vorticity
%option 1: %one Gaussian at the centre
%w = exp(-2*X.^2-Y.^2/20); 

%option 2: %two same gaussians voticies next to each other which can be
%made to collide
%w = exp(-2*(X+2).^2-(Y).^2/20) + exp(-2*(X-2).^2-(Y).^2/20); 

%option 3: %multiple gaussians
% w = exp(-2*(X+6).^2-(Y).^2/10) + exp(-2*(X-6).^2-(Y).^2/10)...
%     +exp(-2*(X).^2-(Y+6).^2/10) + exp(-2*(X).^2-(Y-6).^2/10); 

%option 4: one positive Gaussian and one negative Gaussian
%w = exp(-2*(X+2).^2-(Y+6).^2/20) - exp(-2*(X-2).^2-(Y+6).^2/20); 

%plot the vorticity
% figure(1)
% pcolor(x,y,abs(w));
% shading interp
% colorbar
% colormap jet
% drawnow %to show the plot right away without waitigng for the code to be 
%finished running


            %%% Fourier Transform of the vorticity %%%
wt = fft2(w);
%stacked data to a coloumn to be passed later to ode45
wt2 = reshape(wt,N,1);

                %%%Iterate and integrate over time %%%
%ode45
%solve ODE (it is now ODE not a PDE any more becasue we 
%descretised the the PDE into a punch of ODEs)
[t,wt2sol] = ode45('spc_rhs',tspan,wt2,[],Kderv, KX, KY, nu,Nx,Ny,N);


                       %%% plotting %%%
for j = 1:length(tspan)
    w=ifft2(reshape(wt2sol(j,:),Nx,Ny));
    figure(2)
    %subplot(4,4,j)
    pcolor(x,y,real(w));
    shading interp
    colorbar
    colormap jet
    drawnow
    axis square
    xlabel('x')
    ylabel('y')
end



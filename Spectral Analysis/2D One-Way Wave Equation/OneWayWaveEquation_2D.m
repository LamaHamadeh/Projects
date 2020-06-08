
close all;
clear all;

                %%% Define 2D Discretised spatial %%%            

%Spatial variable on x direction
Lx=2; %domain on x
delta=0.05; %spatial step size
xmin=-Lx; %minimum boundary
xmax=Lx; %maximum boundary 
Nx=(xmax-xmin)/delta; %number of spatial points
x=linspace(xmin,xmax,Nx); %spatial vector

%Spatial variable on y direction
Ly=2; %domain on y
delta=0.05; %spatial step size
ymin=-Ly; %minimum boundary
ymax=Ly; %maximum boundary 
Ny=(ymax-ymin)/delta; %number of spatial points
y=linspace(ymin,ymax,Ny); %spatial vector

%Total matrix size
N = (Nx * Ny);
%--------------------------

                        %%% 2D Initial state %%%

%initial consitions
[X,Y] = meshgrid(x,y); %Nx X Ny grid
%Gaussian
sigma = 0.5;
%U = exp((-X.^2-Y.^2)/sigma^2); %Gaussian at the centre (0,0)
U = exp((-(X-1).^2-(Y-1).^2)/sigma^2); %Gaussian at the edge centered at (1,1)
% %plotting
% figure(1)
% surfl(X,Y,U);
% shading interp
% colorbar, colormap gray
% xlabel('$x$','Interpreter','latex')
% ylabel('$y$','Interpreter','latex')
% zlabel('$U(x,t=0)$','Interpreter','latex')
% title('Initial Gaussian function')
% drawnow
% set(gca,'FontSize',16)
%--------------------------

                    %%% 2D Wave vector disretisation %%%
                   
%x direction
kx = (2*pi/Lx)*[0:(Nx/2-1) (-Nx/2):-1]'; 
kx(1) = 10^(-6);
%y diresction
ky = (2*pi/Ly)*[0:(Ny/2-1) (-Ny/2):-1]'; 
ky(1) = 10^(-6);
%to give kx and ky the sense of direction
[KX,KY] = meshgrid(kx,ky); %N X N grid

%convert to columns so they can pass to ode45
KX = reshape(KX,N,1);
KY = reshape(KY,N,1);
%--------------------------

                    %%% Fast Fourier Transform %%%
                    
%FT+shift of the initial condition
Ut = fftshift(fft2(U));   
%stacked data to a coloumn to be passed later to ode45
Ut = reshape(Ut,N,1);
%--------------------------

                    %%% Time variable %%%

%tspan = [0 0.01 0.05 0.1 0.2 0.5 0.6 0.8 1 1.5 2 2.5 3 3.5 4];
dt = 0.1; %time step
tmin = 0;
tmax = 4;
tspan = [tmin tmax];
%--------------------------

                %%%Iterate and integrate over time %%%
   
 for TimeIteration = 1:100
    t= TimeIteration * dt;               
    %solve
    [Time,Sol] = ode45('FFT_rhs_2D',tspan,Ut,[], KX, KY);
    %inverse of FT
    Sol = ifft2(ifftshift(reshape(Sol(TimeIteration,:),Nx,Nx))); 
    %plotting
    figure(2)
    %surfl
    subplot(1,2,1)
    surfl(abs(Sol));
    colormap jet
    shading interp
    %colorbar
    axis square
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('$|U(x,t)|$','Interpreter','latex')
    xlim ([0 80]);
    ylim ([0 80]);
    zlim ([0, 1]); %This is very important so we can lock the zoom-in behaviour on the z axis
    set(gca,'FontSize',16)
    txt = {['t = ' num2str(t)]};
    text(60,20,0.8,txt,'FontSize',16)
    %imagesc
    subplot(1,2,2)
    imagesc(abs(Sol));
    colormap jet
    shading interp
    colorbar
    axis square
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    xlim ([0 80]);
    ylim ([0 80]);
    set(gca,'FontSize',16)
    
    suptitle('Spectral Solution of the 2D One-Way Wave Equation')
    drawnow;  
    
 end
 %--------------------------

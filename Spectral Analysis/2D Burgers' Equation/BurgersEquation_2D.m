
close all;
clear all;

                %%% Define 1D Discretised spatial %%%            

%Spatial variable on x direction
Lx=3; %domain on x
delta=0.05; %spatial step size
xmin=-Lx; %minimum boundary
xmax=Lx; %maximum boundary 
Nx=(xmax-xmin)/delta; %number of spatial points
x=linspace(xmin,xmax,Nx); %spatial vector

%Spatial variable on y direction
Ly=3; %domain on y
delta=0.05; %spatial step size
ymin=-Ly; %minimum boundary
ymax=Ly; %maximum boundary 
Ny=(ymax-ymin)/delta; %number of spatial points
y=linspace(ymin,ymax,Ny); %spatial vector

%--------------------------
                	     %%% 2D Initial state %%%

%initial consitions
[X,Y] = meshgrid(x,y); %Nx X Ny grid
%IC
sigma = 1;
U = exp(-((X-1).^2+Y.^2)/sigma); %Gaussian
%plotting
% figure(1)
% surfl(X,Y,U);
% shading interp
% colorbar
% colormap jet
% xlabel('$x$','Interpreter','latex')
% ylabel('$y$','Interpreter','latex')
% zlabel('$U(x,t)=\exp\left[-\frac{(x^2+y^2)}{\sigma}\right]$','Interpreter','latex')
% title('Initial Gaussian function')
% set(gca,'FontSize',16)
%--------------------------

                    %%% Time variable %%%

%tspan = [0 0.01 0.05 0.1 0.2 0.5 0.6 0.8 1 1.5 2 2.5 3 3.5 4];
dt = 0.1; %time step
tmin = 0;
tmax = 4;
tspan = [tmin tmax];
% % % %--------------------------
 
                 %%%Iterate and integrate over time %%%

Time = 50; %maximum time to show the dynamics of the shock wave/solution

 for TimeIteration = 1:3:Time
    t= TimeIteration * dt;               
    %solve
    [Time,Sol] = ode45('FFT_rhs_2D',tspan,U,[],Nx,Ny,Lx,Ly);
    Sol = reshape(Sol(TimeIteration,:),Nx,Ny); 
    %plotting
    h = surf(X,Y,abs(Sol));
    colormap jet
    shading interp
    colorbar
    axis square
    xlim ([-Lx Lx]);
    ylim ([-Ly Ly]);
    zlim ([0, 1]); %This is very important so we can lock the zoom-in behaviour on the z axis
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('$|{U(x,t)|}$','Interpreter','latex')
    txt = {['t = ' num2str(t)]};
    text(1,1,1,txt,'FontSize',14)
    axis square
    set(gca,'FontSize',16)
%     direction = [0 1 0];
%     rotate(h,direction,45)
    drawnow;    
 end



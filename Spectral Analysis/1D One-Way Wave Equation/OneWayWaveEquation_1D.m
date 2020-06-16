
close all;
clear all;

                %%% Define 1D Discretised spatial %%%            

%Spatial variable on x direction
L=2; %domain on x
delta=0.05; %spatial step size
xmin=-L; %minimum boundary
xmax=L; %maximum boundary 
N=(xmax-xmin)/delta; %number of spatial points
x=linspace(xmin,xmax,N); %spatial vector
%--------------------------
                	     %%% 1D Initial state %%%

%Gaussian
sigma = 0.5;
%U = exp(-x.^2/sigma^2); %Gaussian
U = sech(5*x); %Hyperbolic
%plotting
figure(1)
plot(x,U,'b','LineWidth',2);
xlabel('$x$','Interpreter','latex')
ylabel('$U(x)$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
title('Initial Gaussian function')
set(gca,'FontSize',16)
%--------------------------
                    %%% 1D Wave vector disretisation %%%
                   
%x direction
k = (2*pi/L)*[0:(N/2-1) (-N/2):-1]'; 
k(1) = 10^(-6);
k = fftshift(k);
%convert to columns so they can pass to ode45
k = reshape(k,N,1);
% %--------------------------
% 
                    %%% Fast Fourier Transform %%%
                    
%FT+shift of the initial condition
Ut = fftshift(fft(U));   
%stacked data to a coloumn to be passed later to ode45
Ut = reshape(Ut,N,1);
% %--------------------------
% 
                    %%% Time variable %%%

%tspan = [0 0.01 0.05 0.1 0.2 0.5 0.6 0.8 1 1.5 2 2.5 3 3.5 4];
dt = 0.1; %time step
tmin = 0;
tmax = 4;
tspan = [tmin tmax];
% %--------------------------
% 
%                 %%%Iterate and integrate over time %%%

DATA = zeros(length(U),100); %initialise the data matrix

 for TimeIteration = 1:100
    t= TimeIteration * dt;               
    %solve
    [Time,Sol] = ode45('FFT_rhs_1D',tspan,Ut,[], k);
    %inverse of FT
    Sol = ifft(ifftshift((Sol(TimeIteration,:)))); 
    %plotting
    figure(2)
    %absolute solution
    subplot(1,2,1)
    plot(x,abs(Sol),'b','LineWidth',2);
    xlabel('$x$','Interpreter','latex')
    ylabel('$|{U(x,t)|}$','Interpreter','latex')
    ylim ([-1.5 1.5])
    xlim ([-3 3])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',16)
%     txt = {['t = ' num2str(t)]};
%     text(1,1,txt,'FontSize',16)
    axis square
    
    subplot(1,2,2)
    waterfall(x,t,abs(Sol))
    xlabel('$x$','Interpreter','latex')
    ylabel('$t$','Interpreter','latex')
    zlabel('$|{U(x,t)|}$','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    colormap jet
    colorbar
    axis square
    ylim ([0 10])
    xlim ([-2 2])

    suptitle({'Absolute Spectral Solution','of the 1D One-Way Wave Equation',['t = ' num2str(t)]})
    set(gca,'FontSize',16)
    
    hold on
    drawnow;
    
    DATA(:,TimeIteration) = abs(Sol); %store the data at each iteration
 end

%plot the 2D data
figure(3)
subplot(1,2,1)
pcolor(DATA)
shading interp
colormap jet
axis square
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
subplot(1,2,2)
surf(DATA)
shading interp
colormap jet
axis square
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
zlabel('$|{U(x,t)|}$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)

%--------------------------

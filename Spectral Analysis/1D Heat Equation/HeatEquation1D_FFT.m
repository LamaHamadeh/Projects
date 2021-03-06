

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
U = exp(-x.^2/sigma^2); %Gaussian
%plotting
figure(1)
plot(x,abs(U));
xlabel('x')
ylabel('U(x)')
title('Initial Gaussian function')
set(gca,'FontSize',16)

% % Square function: here were can see Gibbs phenomenon on the edges
% % indicating the weakness of appling Fourier Transform on such function
% % that have shar edges, i.e., discontinuity
% U = (2/pi)*(atan(exp(-1i*(pi/2).*x))+atan(exp(1i*(pi/2).*x))); 
% % plotting
% figure(1)
% plot(x,U);
% xlabel('x')
% ylabel('U(x)')
% title('Initial Square function')
% set(gca,'FontSize',16)
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
%    
 for TimeIteration = 1:5:500
    t= TimeIteration * dt;               
    %solve
    [Time,Sol] = ode45('FFT_rhs_1D',tspan,Ut,[], k);
    Sol = ifft(ifftshift(Sol(TimeIteration,:))); 
    plot(x,abs(Sol),'b');
    xlabel('$x$','Interpreter','latex')
    ylabel('$|{U(x,t)|}$','Interpreter','latex')
    ylim ([0 1])
    xlim ([-L L])
    title('Spectral Solution of the 1D Diffusion Equation')
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',16)
    txt = {['t = ' num2str(t)]};
    text(1,0.85,txt,'FontSize',16)

    drawnow;
    pause(0.1)
 end
    %--------------------------



%%Still Working on it!

%Solve 1D non-linear Schrodinger equation (Gross Pitaevskii equation)
%using the operator splitting technique and Fast Fourier Transform

clear all;
close all;

% %Spatial variable on x direction
% Lx=5; %domain length
% delta=0.05; %spatial step size
% xmin=-Lx; %minimum boundary
% xmax=Lx; %maximum boundary 
% Nx=(xmax-xmin)/delta; %number of spatial points
% x=linspace(xmin,xmax,Nx); %spatial vector
% %--------------
% 
% %initial non-interacting BEC wavefunction/state
% sigma = 3; %the width of the initial Gaussian
% psi0 = (1/pi)^(1/4).*exp(-x.^2/sigma);
% psi0 = reshape(psi0,Nx,1);
% % %plotting
% % figure;
% % subplot(2,1,1)
% % plot(x,psi0,'k','LineWidth',2)
% % xlabel('$x$','Interpreter','latex')
% % ylabel('$\psi(x,t=0)$','Interpreter','latex')
% % set(gca,'TickLabelInterpreter','latex')
% % set(gca,'FontSize',16)
% % ylim([-1 1.5])
% % xlim([-15 15])
% % hold on
% %--------------
% 
% % 1D Wave vector disretisation                  
% k = (2*pi/Lx)*[0:(Nx/2-1) (-Nx/2):-1];
% k(1) = 10^(-6);
% k = fftshift(k);
% k = reshape(k,Nx,1);
% %--------------
% 
% %FFT of IC
% % psi0t = fftshift(fft(psi0));
% % %plotting
% % subplot(2,1,2)
% % plot(k,abs(psi0t)/max(psi0t),'k','LineWidth',2)
% % xlabel('$k$','Interpreter','latex')
% % ylabel('$\hat{\psi}(k,t=0)$','Interpreter','latex')
% % set(gca,'TickLabelInterpreter','latex')
% % set(gca,'FontSize',16)
% % ylim([-1 1.5])
% % xlim([-15 15])
% %--------------
% 
% %Time variable
% dt = 0.1; %time step
% tmin = 0;
% tmax = 4;
% tspan = [tmin tmax];
% %--------------
% 
% %solve
% 
% % initialise space and time data 
% psidata = psi0; 
% tdata = 0; 
% 
% Time = 20; %maximum time
% 
% for i =1:Time
%     t = i * dt;
%     psi0t = fftshift(fft(psi0));
%     psi1 = ifft(ifftshift(psi0t.*exp(-1i.*k.^2*t/2))); %the solution of the dispersion term
%     psi2 = psi1.*exp(1i.*(abs(psi1).^2)*t);
%     psi0 = psi2;
%     
%     u = abs(psi0);
%     psidata = [psidata u]; %space matrix (space x time)
%     tdata = [tdata t]; %time row matrix (1 x time)
%         
%     plot(x,psidata(:,i))
%     xlabel('$x$','Interpreter','latex')
%     ylabel('$|\psi(x,t)|$','Interpreter','latex')
%     set(gca,'TickLabelInterpreter','latex')
%     set(gca,'FontSize',16)
%     xlim([-7 7])
%     ylim([-1 1.5])
%     
%     drawnow
%     pause(0.2)
%     
% end
% 
% 
% % %plotting
% % figure
% % waterfall(x,tdata,psidata')
% % colormap jet
% % % surf(x,tdata,psidata')
% % % shading interp
% % xlabel('$x$','Interpreter','latex')
% % ylabel('$t$','Interpreter','latex')
% % zlabel('$|\psi(x,t)|$','Interpreter','latex')
% % set(gca,'TickLabelInterpreter','latex')
% % set(gca,'FontSize',16)
% % grid off
% % axis square


%Non-linear Schodinger Equation
%-------------------------------

%Leap Frog Scheme (2,2)
%----------------------

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

CFL = dt/dx^2;

%initial condition
%this is a two-step procedure therefore we nedd two initial conditions.
u0 = exp(-x.^2).'; %the transpose is to set the values as a column
u1 = exp(-(x).^2).';
usol(:,1) = u0; 
usol(:,2) = u1;
%plot(x,u0);
 
%u_xx matrix
e1 = ones(n,1);
A = spdiags([e1 -2*e1 e1], [-1 0 1], n, n);
A(1,n) = 1;
A(n,1) = 1;
 
%set up the iteration process 
for j = 1:length(t)-2
    u2 = u0 -(1i*CFL*A*u1-1i*2*dt*(conj(u1).*u1).*u1);
    u0 = u1; %my past is my current
    u1 = u2; %my current is my future
    usol(:,j+2) = u2;
end
% 
%plot the one way wave
figure(5)
surf(x,t,abs(usol'));
% map = [ 0 0 0];
colormap jet;
xlabel('Distance')
ylabel('Time')
zlabel('Wave Solution')



%ideas
%1- solve 1D time depdnent schrdoinger for non interacting BEC gas in a
%harmonic oscillator trap (investigate egienstates and eigenvalues and the evolution of the initial state)
%2- solve 1D time depdnent schrdoinger for non interacting BEC gas in a
%double-well trap (investigate egienstates and eigenvalues and the
%evolution of the initial state). Here it is important to include the
%analysis of the quantum tunelling
%introduce a time-depdnent trap: from harmonic to double well
%introduce the interacting gas and the gross pitaevskii equation. 

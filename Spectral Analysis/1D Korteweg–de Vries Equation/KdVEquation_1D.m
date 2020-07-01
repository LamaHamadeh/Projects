%%%Still working on it%%%

close all;
clear all;

                %%% Define 1D Discretised spatial %%%            

%Spatial variable on x direction
L=4; %domain on x
delta=0.05; %spatial step size
xmin=-L; %minimum boundary
xmax=L; %maximum boundary 
N=(xmax-xmin)/delta; %number of spatial points
x=linspace(xmin,xmax,N); %spatial vector
%--------------------------
                	     %%% 1D Initial state %%%

%IC
% sigma = 0.5;
% U = exp(-x.^2/sigma^2); %Gaussian
% U = cos(pi.*x);
% A = 25; 
% B = 16;
% U = A^2*sech(.5*(A*(x+2))).^2 + B^2*sech(.5*(B*(x+1))).^2;
sigma1 = 2; sigma2 = 4;
c1     = 2;     c2 = 1;
h1     = 1;     h2 = 0.3;
U = 0.*x;
U = h1*sech(sigma2*(x+c1)).^2+h2*sech(sigma1*(x-c2)).^2;
% %plotting
% figure(1)
% plot(x,U,'b','LineWidth',2);
% xlabel('$x$','Interpreter','latex')
% ylabel('$U(x,t=0)$','Interpreter','latex')
% xlim([-L-2 L+2])
% ylim([-0.5 1.2])
% set(gca,'TickLabelInterpreter','latex')
% title('Initial State')
% set(gca,'FontSize',16)
% drawnow;
%--------------------------

                    %%% Time variable %%%

%tspan = [0 0.01 0.05 0.1 0.2 0.5 0.6 0.8 1 1.5 2 2.5 3 3.5 4];
dt = 0.01; %time step
tmin = 0;
tmax = 4;
tspan = [tmin tmax];
% %--------------------------

%                 %%%Iterate and integrate over time %%%

% DATA = zeros(length(U),100); %initialise the data matrix
Time = 100; 

 for TimeIteration = 1:2:Time
    t= TimeIteration * dt;               
    %solve
    [Time,Sol] = ode113('FFT_rhs_1D',tspan,U,[],N,L); %deals with the
%   stiffness by using implicit schemes/ predictor - corrector (or
%   ode15s - Gear methods), or ode23.
    %inverse of FT
    Sol = Sol(TimeIteration,:);
    %plotting
    figure(2)
    %absolute solution
    subplot(1,2,1)
    plot(x,abs(Sol),'b','LineWidth',2);
    xlabel('$x$','Interpreter','latex')
    ylabel('$|{U(x,t)|}$','Interpreter','latex')
    xlim([-L-2 L+2])
    ylim([-0.5 1.2])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',16)
    txt = {['t = ' num2str(t)]};
    text(1,1,txt,'FontSize',16)
    axis square
    
%     subplot(1,2,2)
%     waterfall(x,t,abs(Sol))
%     xlabel('$x$','Interpreter','latex')
%     ylabel('$t$','Interpreter','latex')
%     zlabel('$|{U(x,t)|}$','Interpreter','latex')
%     set(gca,'TickLabelInterpreter','latex')
%     colormap jet
%     colorbar
%     axis square
%     ylim ([0 10])
%     xlim ([-L L])
% 
%     suptitle({'Absolute Spectral Solution','of the 1D KdV Equation',['t = ' num2str(t)]})
%     set(gca,'FontSize',16)
    
    drawnow
%     DATA(:,TimeIteration) = abs(Sol); %store the data at each iteration
 end
% % 
% % %plot the 2D data
% % figure(3)
% % subplot(1,2,1)
% % pcolor(DATA)
% % shading interp
% % colormap jet
% % axis square
% % xlabel('$t$','Interpreter','latex')
% % ylabel('$x$','Interpreter','latex')
% % set(gca,'TickLabelInterpreter','latex')
% % set(gca,'FontSize',16)
% % subplot(1,2,2)
% % surf(DATA)
% % shading interp
% % colormap jet
% % axis square
% % xlabel('$t$','Interpreter','latex')
% % ylabel('$x$','Interpreter','latex')
% % zlabel('$|{U(x,t)|}$','Interpreter','latex')
% % set(gca,'TickLabelInterpreter','latex')
% % set(gca,'FontSize',16)
% 
% 
% %Explanation:
% % Korteweg-de Vries equation (KDVE) is a nonlinear partial differential 
% %equation, which is initially introduced to describe the lossless 
% %propagation of shallow water waves. It represents the long time 
% %evolution of wave phenomena, in which the effect of the nonlinear 
% %terms UUx is counterbalanced by the dispersion Uxxx.Thus it has been found 
% %to model many wave phenomena such as waves in enharmonic crystals, bubble 
% %liquid mixtures, ion acoustic wave and magnetohydrodynamic waves in a warm 
% %plasma as well as shallow water waves.
% 
% % Among many remarkable properties (all stemming from the fact that it is 
% %"integrable" in some sense), the nonlinear term (the second term on the 
% %right hand side) somehow exactly compensates for the dispersive term 
% %(the first term on the right hand side). This allows it to admit "soliton"
% %solutions.
% 
% % Furthermore, higher solitons travel faster. This means a tall soliton 
% %will catch up with a short one, and they will undergo some strange 
% %interaction. Amazingly, despite the nonlinearity, a long time afterwards 
% %they seem to recover their original forms (apart from a phase shift). 
% %In other words, the tall one seems to have overtaken the short one!
% % 
% % This preservation of structure suggests that there are many quantities 
% %conserved under evolution, and in fact, there are infinitely many.
% 
% %--------------------------
% % 
% % %Spectral Methods in Matlab By Lioyd N. Trefethen
% % 
% % % p27.m - Solve KdV eq. u_t + uu_x + u_xxx = 0 on [-pi,pi] by
% % % FFT with integrating factor v = exp(-ik^3t)*u-hat.
% % % Set up grid and two-soliton initial data:
% % N = 256; dt = .4/N^2; x = (2*pi/N)*(-N/2:N/2-1)';
% % A = 25; B = 16; clf, drawnow
% % u = 3*A^2*sech(.5*(A*(x+2))).^2 + 3*B^2*sech(.5*(B*(x+1))).^2;
% % v = fft(u); k = [0:N/2-1 0 -N/2+1:-1]'; ik3 = 1i*k.^3;
% % % Solve PDE and plot results:
% % tmax = 0.006; nplt = floor((tmax/25)/dt); nmax = round(tmax/dt);
% % udata = u; tdata = 0; h = waitbar(0,'please wait...');
% % for n = 1:nmax
% % t = n*dt; g = -.5i*dt*k;
% % E = exp(dt*ik3/2); E2 = E.^2;
% % a = g.*fft(real( ifft( v ) ).^2);
% % b = g.*fft(real( ifft(E.*(v+a/2)) ).^2); % 4th-order
% % c = g.*fft(real( ifft(E.*v + b/2) ).^2); % Runge-Kutta
% % d = g.*fft(real( ifft(E2.*v+E.*c) ).^2);
% % v = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6;
% % if mod(n,nplt) == 0
% % u = real(ifft(v)); waitbar(n/nmax)
% % udata = [udata u]; tdata = [tdata t];
% % end
% % end
% % waterfall(x,tdata,udata'), colormap([0 0 0]), view(-20,25)
% % xlabel x, ylabel t, axis([-pi pi 0 tmax 0 2000]), grid off
% % set(gca,'ztick',[0 2000]), close(h), pbaspect([1 1 .13])
% % %--------------------------




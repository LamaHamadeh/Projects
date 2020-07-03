
%Code adopted from "Spectral Methods in MATLAB" By Lloyd N. Trefethen 
% code: p27.m, page: 110, 111, 112, 113
%with some modification and added comments

% Solve 1D KdV eq. u_t + uu_x + u_xxx = 0 on [-pi,pi] by
% FFT with integrating factor v = exp(-ik^3t)*u-hat.
% Or: spectral filtering to remove the stiffness resulting
%from the high order linear term.

% Set up space grid
N = 256;
x = (2*pi/N)*(-N/2:N/2-1)';
nmax = round(tmax/dt); %number of space points

% Two solitons as initial condition
A = 25; 
B = 16; 
u = x.*0;
u = 3*A^2*sech(.5*(A*(x+2))).^2 + 3*B^2*sech(.5*(B*(x+1))).^2;

% Fourier Transform of IC
v = fft(u); 
U = v; 

% Define wave number variable
k = [0:N/2-1 0 -N/2+1:-1]'; 
ik3 = 1i*k.^3;

% Define time variable
dt = 0.4/N^2; %time step
tmax = 0.006; %maximum time 
nplt = floor((tmax/25)/dt); %number of time points

% initialise space and time data 
udata = u; 
tdata = 0; 

% Solve PDE and plot results:
for n = 1:nmax
   
    % Iterate each time step
    t = n*dt;
    
    % Define common factor appearing for each RK4 terms
    g = -.5i*dt*k; 
    
    % Define exponentials
    Etp = exp(t*ik3); %posiitve
    Etm = exp(-t*ik3); %negative
    %positive time step
    Ep = exp(dt*ik3/2); %half
    Ep2 = Ep.^2; %full
    %negative time step
    Em = exp(-dt*ik3/2); %half
    Em2 = Em.^2; %full

    % 4th-order Runge-Kutta Method
    % F1 = dt f(yn,tn)
    F1 = g.* Etm.* fft(( ifft(Etp.* U ) ).^2); %F1
    % F2 = dt f(yn+F1/2,tn+dt/2)
    F2 = g.* Etm .* Em .* fft(( ifft(Etp.*Ep.*(U+F1/2)) ).^2); %F2
    % F3 = dt f(yn+F2/2,tn+dt/2)
    F3 = g.* Etm .* Em .* fft(( ifft(Etp.*Ep.*(U+F2/2)) ).^2); %F3
    % F4 = dt f(yn+F3,tn+dt)
    F4 = g.* Etm .* Em2 .* fft(( ifft(Etp.*Ep2.*(U+F3)) ).^2); %F4
    % The approximated future function after one time step: yn+1 = yn + (1/6)*[F1+2(F2+F3)+F4]
    U = U + (F1+2*(F2+F3)+F4)/6;

%     if mod(n,nplt) == 0
        %To get back to our actual solution: uhat = Uhat*exp(i*K^3*t)
        u = abs(ifft(U.*Etp));
        %store the space and time data at each iteration
        udata = [udata u]; %space matrix (space x time)
        tdata = [tdata t]; %time row matrix (1 x time)
%     end
   
end

%plotting
figure
% subplot(1,2,1)
waterfall(x,tdata,udata')
colormap(jet)
view(-20,25)
xlabel('$x$','Interpreter','latex')
ylabel('$t$','Interpreter','latex')
zlabel('$|u(x,t)|$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
axis([-pi pi 0 tmax 0 2000])
grid off
axis square
set(gca,'ztick',[0 2000])
pbaspect([1 1 .13])

% subplot(1,2,2)
% imagesc(tdata,x,udata)
% colormap(jet)
% xlabel('$t$','Interpreter','latex')
% ylabel('$x$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'FontSize',16)
% axis square





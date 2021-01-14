
close all;
clear all;

%Nonlinear partial differential Schrodinger Equation
%iu_t + 1/2 u_xx + |u|^2 u = 0
%with u(0,x) = 2sech(x).


%Spectral solution of NLS

%space
L = 30;
n = 1024;
x2 = linspace(-L/2, L/2, n+1);
x = x2(1:n); %periodic so cut off endpoint.
k = (2*pi/L) * [0:n/2-1 -n/2:-1].'; %wave number 

%time
slices = 40;
t = linspace(0, 2*pi, slices+1);
dt = t(2) - t(1);

%initial conditions
u0 = 2*sech(x).';

%solve in frequency domain
ut = fft(u0);
[t, utsol] = ode45('nls_rhs', t, ut, [], k);

%bring back to space
usol = zeros(length(t),n);
for j=1:length(t)
    usol(j,:) = ifft(utsol(j,:));
end

%The solution matrix
X = usol.';
[md,nd] = size(X);

%visualise NLSE solution
figure(1)
waterfall(x,t,abs(usol))
colormap([0,0,0])
set(gca, 'Ylim', [0,2*pi],'Xlim',[-L/2 L/2])
zlabel('$\left|u\right|$','interpreter','latex')
xlabel('$x$','interpreter','latex')
ylabel('$t$','interpreter','latex')
title('Numerical Solution of the NLSE')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)


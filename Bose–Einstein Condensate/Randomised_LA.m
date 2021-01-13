
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
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';

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

X = usol.';
[md,nd] = size(X);

figure(1)
waterfall(x,t,abs(usol))
colormap([0,0,0])
set(gca, 'Ylim', [0,2*pi],'Xlim',[-L/2 L/2])
title('|u|')

%apply singualr value decomposition
[u,s,v] = svd(X, 'econ');
%show the eigenvalues
figure(2)
plot(diag(s)/sum(diag(s)),'ko','Linewidth',2); %gives the percetages of each mode
figure(3)
semilogy(diag(s)/sum(diag(s)),'ko','Linewidth',2)
%show the spatial dominant modes
figure(4)
plot(x,u(:,1),x,u(:,2),x,u(:,3),'Linewidth',2)
%show the temporal dominant modes
figure(5)
plot(t,v(:,1),t,v(:,2),t,v(:,3),'Linewidth',2)

%------------------------

%randomised LA/projections
%STAGE A: Take the data matrix, do random projection to contruct a set of
%orthonormal vectors that span the column space of the data matrix 
%STEPS 1: random projection
[M,N] = size(X);
%number of random projections I pick (randomised sub-sampling)
K = 10; 
Omega=randn(N,K);
Y = X*Omega;
size(Y) 
%STEP 2: DQ Decompositon
[Q,R] = qr(Y,0);
size(Q)

%STAGE B
B = (Q')*X;
%randomised approximations
[U,S,V] = svd(B,'econ');
%Go back to the full decompositon as an approximation
uapprox = Q*U;
%show the spatial dominant modes: actual and approximations
figure(6)
subplot(3,1,1)
plot(x,u(:,1),'k',x,uapprox(:,1),'k.','Linewidth',2)
subplot(3,1,2)
plot(x,u(:,2),'b',x,uapprox(:,2),'b.','Linewidth',2)
subplot(3,1,3)
plot(x,u(:,3),'r',x,uapprox(:,3),'r.','Linewidth',2)

%------------------------




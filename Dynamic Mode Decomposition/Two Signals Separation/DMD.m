clear all;
close all;


%% spatiotemporal variables
%space
xi = linspace(-10,10,400);

%time
t = linspace(0,4*pi,200);
dt = t(2)-t(1); %time step

%spatio-temporal grid
[Xgrid,T] = meshgrid(xi,t);

%% create two spatio-temporal patterns, numerical data
f1 = sech(Xgrid+3).*(1*exp(1i*2.3*T));
f2 = (sech(Xgrid).*tanh(Xgrid)).*(2*exp(1i*2.8*T));
f = f1+f2;
%first function
figure(1)
subplot(2,2,1)
surfl(Xgrid,T,real(f1));
shading interp
colormap (gray)
%second function
subplot(2,2,2)
surfl(Xgrid,T,real(f2));
shading interp
colormap (gray)
%sum of both function
subplot(2,2,3)
surfl(Xgrid,T,real(f));
shading interp
colormap (gray)

%% svd, low rank structure
[u,s,v] = svd(f.');
%eigenvalues of svd
figure(2)
plot(diag(s)/(sum(diag(s))),'ro')

%plot 
figure(3)
subplot(2,1,1)
plot(real(u(:,1:2)),'Linewidth',2)%modes

%plot
subplot(2,1,2)
plot(real(v(:,1:2))) %what modes behave as a function of time

%% DMD
X = f.'; %rows are spatial data and columns are time snapshots.

X1 = X(:,1:end-1);
X2 = X(:,2:end);

%rank trauncation (2 modes)
r = 2;
[U,S,V] = svd(X1,'econ');

Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);

%least square matrix
Atilde = Ur'*X2*Vr/Sr;
[W,D] = eig(Atilde); %eigenvectors and eigenvalues

Phi = X2*Vr/Sr*W; %DMD modes

lambda = diag(D);
omega = log(lambda)/dt;

figure(3)
subplot(2,1,1), hold on
plot(real(Phi),'Linewidth',2)


%reconstruction of modes/function in time
X1 = X(:,1); %t=0
b = Phi\X1;

%future preiction
t2 = linspace(0,20*pi,200);

time_dynamics = zeros(r,length(t));
for iter = 1:length(t)
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end

X_DMD = Phi * time_dynamics;
figure(1)
subplot(2,2,4)
surfl(Xgrid,T,real(X_DMD).')
shading interp
colormap gray

%check the eignevalues that all of them imaginary
figure(4)
plot(real(omega),imag(omega),'r*','Linewidth',2)
xlabel('Real')
ylabel('Imaginary')
xlim([0 0.2])
ylim([2 3])







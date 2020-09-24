clear all;
close all;

%%%Define numerical discretisation variables
%%%-----------------------------------------
%space
L=30;
n=512;
x2=linspace(-L/2,L/2,n+1);
x=x2(1:n);
%wave number (since we will use spectral method to solve our PDE)
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
%time
slices=20;
t=linspace(0,pi,slices+1);
dt=t(2)-t(1);

%%%%NLS equation solve (creation of data)
%%%-----------------------------------------
%initial conditions
N=2;
u=N*(sech(x)).';
%Fourier transform of initial condition
ut=fft(u);
%Runge Kutta solver
[t,utsol]=ode45('NLSE_solution_rhs',t,ut,[],k);
%inverse of Fourier transform and back to time domain on every time step
for j=1:length(t)
    usol(j,:)=ifft(utsol(j,:));  % bring back to space
end
%define the spatio-temporal data
X = usol.';  

%%%%DMD
%%%-----
%define the DMD matrices
X1 = X(:,1:end-1);
X2 = X(:,2:end);
%SVD
[U1,Sigma1,V1] = svd(X1, 'econ');
%define the DMD truncation rank
rDMD=10;
%define the truncated matrices
UDMD=U1(:,1:rDMD);
SigmaDMD=Sigma1(1:rDMD,1:rDMD);
VDMD=V1(:,1:rDMD);
%define the linear operator
Atilde_DMD = UDMD'*X2*VDMD/SigmaDMD;
%eigen decomposition of the linear operator
[WDMD,DDMD] = eig(Atilde_DMD);
%DMD eigenvector
Phi=X2*VDMD/SigmaDMD*WDMD;
%DMD eigenvalues (dynamics frequencies)
lambda1=diag(DDMD);
OmegaDMD=log(lambda1)/dt;
% pseudo-inverse initial conditions
y01 = Phi\u;  
% DMD reconstruction for every time point
u_modes = zeros(rDMD,length(t));  
for iter = 1:length(t)
    u_modes(:,iter) =(y01.*exp(OmegaDMD*(t(iter))));
end
% DMD resconstruction with all modes
u_dmd = Phi*u_modes;   

%error estimation
for j=1:length(t)
    error_DMD(j)=norm(u_dmd(:,j)-usol(j,:).');
end

%%%%Koopman (u.*|u|^2 as an observable)
%%%------------------------------------
%define the DMD matrices
Y1=[X1; (X1.*abs(X1).^2)];
Y2=[X2; (X2.*abs(X2).^2)];
%SVD
[U2,Sigma2,V2] = svd(Y1, 'econ');
%define the DMD truncation rank
rKoop=10; 
%define the truncated matrices
UKoop=U2(:,1:rKoop); 
SigmaKoop=Sigma2(1:rKoop,1:rKoop); 
VKoop=V2(:,1:rKoop);
%define the linear operator
Atilde_Koop = UKoop'*Y2*VKoop/SigmaKoop;
%eigen decomposition of the linear operator
[WKoop,DKoop] = eig(Atilde_Koop);
%DMD eigenvector
Phi2=Y2*VKoop/SigmaKoop*WKoop;
%DMD eigenvalues (dynamics frequencies)
lambda2=diag(DKoop);
OmegaKoop=log(lambda2)/dt;
%define a vector that has the same size of Phi2
u2=[u; (u.*abs(u).^2)];
% pseudo-inverse initial conditions
y0 = Phi2\u2;  
% DMD reconstruction for every time point
u_modes = zeros(rKoop,length(t));  % DMD reconstruction for every time point
for iter = 1:length(t)
    u_modes(:,iter) =(y0.*exp(OmegaKoop*(t(iter))));
end
% DMD resconstruction with all modes
u_dmd2 = Phi2*u_modes;   
%define the resulting Koopman modes within the limits of our dyanmics
u_Koop = u_dmd2(1:n,:);
Phi_Koop = Phi2(1:n,:);

%error estimation
for j=1:length(t)
    error_Koop1(j)=norm(u_Koop(:,j)-usol(j,:).');
end


%%%%Koopman (|u|^2 as an observable)
%%%------------------------------------
%define the DMD matrices
W1=[X1; (abs(X1).^2)];
W2=[X2; (abs(X2).^2)];
%SVD
[U3,Sigma3,V3] = svd(W1, 'econ');
%define the DMD truncation rank
rKoop2=10; 
%define the truncated matrices
UKoop2=U3(:,1:rKoop2); 
SigmaKoop2=Sigma3(1:rKoop2,1:rKoop2); 
VKoop2=V3(:,1:rKoop2);
%define the linear operator
Atilde_Koop2 = UKoop2'*W2*VKoop2/SigmaKoop2;
%eigen decomposition of the linear operator
[WKoop2,DKoop2] = eig(Atilde_Koop2);
%DMD eigenvector
Phi3=W2*VKoop2/SigmaKoop2*WKoop2;
%DMD eigenvalues (dynamics frequencies)
lambda3=diag(DKoop2);
OmegaKoop2=log(lambda3)/dt;
%define a vector that has the same size of Phi2
u2=[u; (u.*abs(u).^2)];
% pseudo-inverse initial conditions
y0 = Phi3\u2;  
% DMD reconstruction for every time point
u_modes = zeros(rKoop2,length(t));  % DMD reconstruction for every time point
for iter = 1:length(t)
    u_modes(:,iter) =(y0.*exp(OmegaKoop2*(t(iter))));
end
% DMD resconstruction with all modes
u_dmd2 = Phi3*u_modes;   
%define the resulting Koopman modes within the limits of our dyanmics
u_Koop2 = u_dmd2(1:n,:);
Phi_Koop2 = Phi2(1:n,:);

%error estimation
for j=1:length(t)
    error_Koop2(j)=norm(u_Koop2(:,j)-usol(j,:).');
end

%%%%Visualisation
%%%--------------
%numerical solution (PDE)
subplot(3,4,1), waterfall(x,t,abs(usol)), colormap([0 0 0])
set(gca,'Ylim',[0 pi],'Ytick',[0 3 6],'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter','latex'), zlabel('$\left|u\right|^2$','interpreter','latex')
axis square, title('PDE','interpreter','latex'),set(gca,'fontsize',14), set(gca,'TickLabelInterpreter','latex')
%DMD
subplot(3,4,2), waterfall(x,t,abs(u_dmd).'), colormap([0 0 0])
set(gca,'Ylim',[0 pi],'Ytick',[0 3 6],'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter','latex'), zlabel('$\left|u\right|^2$','interpreter','latex')
axis square, title('DMD','interpreter','latex'),set(gca,'fontsize',14), set(gca,'TickLabelInterpreter','latex')
%Koop (u.|u|^2)
subplot(3,4,3), waterfall(x,t,abs(u_Koop).'), colormap([0 0 0])
set(gca,'Ylim',[0 pi],'Ytick',[0 3 6],'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter','latex'), zlabel('$\left|u\right|^2$','interpreter','latex')
axis square, title('Koopman $\left(u.|u|^2\right)$','interpreter','latex')
set(gca,'fontsize',14), set(gca,'TickLabelInterpreter','latex')
%Koop (|u|^2)
subplot(3,4,4), waterfall(x,t,abs(u_Koop2).'), colormap([0 0 0])
set(gca,'Ylim',[0 pi],'Ytick',[0 3 6],'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 10],'Ztick',[0 5 10])
xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter','latex'), zlabel('$\left|u\right|^2$','interpreter','latex')
axis square, title('Koopman $\left(|u|^2\right)$','interpreter','latex')
set(gca,'fontsize',14), set(gca,'TickLabelInterpreter','latex')

%eigenvalues
%DMD
subplot(3,4,6)
plot(OmegaDMD,'ko','Linewidth',2)
grid on
axis square
axis([-5 1 -20 20])
xlabel('Real')
ylabel('Imaginary')
set(gca,'fontsize',14)
set(gca,'TickLabelInterpreter','latex')
%Koop (u.|u|^2)
subplot(3,4,7)
plot(OmegaKoop,'bo','Linewidth',2)
grid on
axis square
axis([-5 1 -20 20])
xlabel('Real')
ylabel('Imaginary')
set(gca,'fontsize',14)
set(gca,'TickLabelInterpreter','latex')
%Koop (|u|^2)
subplot(3,4,8)
plot(OmegaKoop2,'ro','Linewidth',2)
grid on
axis square
axis([-5 1 -20 20])
xlabel('Real')
ylabel('Imaginary')
set(gca,'fontsize',14)
set(gca,'TickLabelInterpreter','latex')

%Error estimaton
subplot(3,4,[9,12])
semilogy(t,error_DMD,'k',t,error_Koop1,'b',t,error_Koop2,'r','LineWidth',2)
set(gca,'Ylim',[10^(-6) 10^2],'Ytick',[10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
xlim([0 max(t)])
xlabel('Time')
ylabel('Error')
set(gca,'fontsize',14)
set(gca,'TickLabelInterpreter','latex')
leg1 = legend('DMD','Koopman $\left(u.|u|^2\right)$','Koopman $\left(|u|^2\right)$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',14);


    

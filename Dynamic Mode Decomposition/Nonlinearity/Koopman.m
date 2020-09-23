clear all;
close all;

%NLS SOLVE - CREATION OF DATA
% space
L=30;
n=512;
x2=linspace(-L/2,L/2,n+1);
x=x2(1:n);
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
% time
slices=20;
t=linspace(0,pi,slices+1);
dt=t(2)-t(1);

% initial conditions
N=2;
u=N*(sech(x)).';
ut=fft(u);
[t,utsol]=ode45('NLSE_solution_rhs',t,ut,[],k);
for j=1:length(t)
    usol(j,:)=ifft(utsol(j,:));  % bring back to space
end
X = usol.';  % here is the data

%--------------

%DMD
X1 = X(:,1:end-1);
X2 = X(:,2:end);

[U1,Sigma1,V1] = svd(X1, 'econ');
rDMD=10;
UDMD=U1(:,1:rDMD);
SigmaDMD=Sigma1(1:rDMD,1:rDMD);
VDMD=V1(:,1:rDMD);
Atilde_DMD = UDMD'*X2*VDMD/SigmaDMD;
[WDMD,DDMD] = eig(Atilde_DMD);
Phi=X2*VDMD/SigmaDMD*WDMD;

%eigenvalues
lambda1=diag(DDMD);
OmegaDMD=log(lambda1)/dt;

y01 = Phi\u;  % pseudo-inverse initial conditions

u_modes = zeros(rDMD,length(t));  % DMD reconstruction for every time point
for iter = 1:length(t)
    u_modes(:,iter) =(y01.*exp(OmegaDMD*(t(iter))));
end
u_dmd = Phi*u_modes;   % DMD resconstruction with all modes

%error estimation
for j=1:length(t)
    error2(j)=norm(u_dmd(:,j)-usol(j,:).');
end

%--------------

%Koopman1
Y1=[X1; (X1.*abs(X1).^2)];
Y2=[X2; (X2.*abs(X2).^2)];
[U2,Sigma2,V2] = svd(Y1, 'econ');
rKoop=10; 
UKoop=U2(:,1:rKoop); 
SigmaKoop=Sigma2(1:rKoop,1:rKoop); 
VKoop=V2(:,1:rKoop);
Atilde_Koop = UKoop'*Y2*VKoop/SigmaKoop;
[WKoop,DKoop] = eig(Atilde_Koop);
Phi2=Y2*VKoop/SigmaKoop*WKoop;

lambda2=diag(DKoop);
OmegaKoop=log(lambda2)/dt;

u2=[u; (u.*abs(u).^2)];
y0 = Phi2\u2;  % pseudo-inverse initial conditions

u_modes = zeros(rKoop,length(t));  % DMD reconstruction for every time point
for iter = 1:length(t)
    u_modes(:,iter) =(y0.*exp(OmegaKoop*(t(iter))));
end
u_dmd2 = Phi2*u_modes;   % DMD resconstruction with all modes

u_Koop = u_dmd2(1:n,:);
Phi_Koop = Phi2(1:n,:);

%error estimation
for j=1:length(t)
    error3(j)=norm(u_Koop(:,j)-usol(j,:).');
end

%--------------

%Koopman2
W1=[X1; (abs(X1).^2)];
W2=[X2; (abs(X2).^2)];
[U3,Sigma3,V3] = svd(W1, 'econ');
rKoop2=10; 
UKoop2=U3(:,1:rKoop2); 
SigmaKoop2=Sigma3(1:rKoop2,1:rKoop2); 
VKoop2=V3(:,1:rKoop2);
Atilde_Koop2 = UKoop2'*W2*VKoop2/SigmaKoop2;
[WKoop2,DKoop2] = eig(Atilde_Koop2);
Phi3=W2*VKoop2/SigmaKoop2*WKoop2;

lambda3=diag(DKoop2);
OmegaKoop2=log(lambda3)/dt;

u2=[u; (u.*abs(u).^2)];
y0 = Phi3\u2;  % pseudo-inverse initial conditions

u_modes = zeros(rKoop2,length(t));  % DMD reconstruction for every time point
for iter = 1:length(t)
    u_modes(:,iter) =(y0.*exp(OmegaKoop2*(t(iter))));
end
u_dmd2 = Phi3*u_modes;   % DMD resconstruction with all modes

u_Koop2 = u_dmd2(1:n,:);
Phi_Koop2 = Phi2(1:n,:);

%error estimation
for j=1:length(t)
    error4(j)=norm(u_Koop2(:,j)-usol(j,:).');
end

%--------------

%visualisation
%numerical solution
subplot(3,4,1), waterfall(x,t,abs(usol)), colormap([0 0 0])
set(gca,'Ylim',[0 pi],'Ytick',[0 3 6],'Fontsize',14,'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('\xi','Fontsize',14), ylabel('time','Fontsize',14)
axis square, title('PDE')
%DMD
subplot(3,4,2), waterfall(x,t,abs(u_dmd).'), colormap([0 0 0])
set(gca,'Ylim',[0 pi],'Ytick',[0 3 6],'Fontsize',14,'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('\xi','Fontsize',14), ylabel('time','Fontsize',14)
axis square, title('DMD')
%Koop
subplot(3,4,3), waterfall(x,t,abs(u_Koop).'), colormap([0 0 0])
set(gca,'Ylim',[0 pi],'Ytick',[0 3 6],'Fontsize',14,'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('\xi','Fontsize',14), ylabel('time','Fontsize',14)
axis square, title('Koopman (u.|u|^2)')
%Koop
subplot(3,4,4), waterfall(x,t,abs(u_Koop2).'), colormap([0 0 0])
set(gca,'Ylim',[0 pi],'Ytick',[0 3 6],'Fontsize',14,'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 10],'Ztick',[0 5 10])
xlabel('\xi','Fontsize',14), ylabel('time','Fontsize',14)
axis square, title('Koopman (|u|^2)')

%eigenvalues
%DMD
subplot(3,4,6)
plot(OmegaDMD,'ko','Linewidth',2), grid on, axis square
%Koop
subplot(3,4,7)
plot(OmegaKoop,'ko','Linewidth',2), grid on, axis square
%Koop
subplot(3,4,8)
plot(OmegaKoop2,'ko','Linewidth',2), grid on, axis square

%error estimaton
subplot(3,4,[9,12])
semilogy(t,error2,'k:',t,error3,'k-.',t,error4,'ks','LineWidth',2)
%set(gca,'Ylim',[10^(-6) 10^2],'Ytick',[10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
xlabel('Time')
ylabel('Error')
legend('DMD', 'Koopman','Koopman')

    

clear all;
close all;


%% Define time and space discretisations
%space
xi = linspace(-10,10,400);
%time
t = linspace(0,4*pi,200);
dt = t(2)-t(1); %time step
%spatio-temporal grid
[Xgrid,T] = meshgrid(xi,t);

%% create two spatio-temporal patterns, numerical data
%first signal
f1 = sech(Xgrid+3).*(1*exp(1i*2.3*T));
%second signal
f2 = (sech(Xgrid).*tanh(Xgrid)).*(2*exp(1i*2.8*T));

%visualise spatial modes
figure;
plot(xi,sech(xi+3),'k','LineWidth',2)
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',16)

figure;
plot(xi,-2.*(sech(xi).*tanh(xi)),'k','LineWidth',2)
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',16)


%combine signals and make data matrix
f = f1+f2;
X = f.'; %Data matrix (%rows are spatial data and columns are time snapshots.)

%visualise f1, f2 and f
%f1
figure,
subplot(2,2,1)
surfl(real(f1));
shading interp
colormap (gray)
view(-20,60)
set(gca , 'YTick', numel(t)/4 * (0:4)),
set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca , 'XTick', linspace(1,numel(xi),3)),
set(gca , 'Xticklabel',{'-10', '0', '10'});
title('$f_1(x,t)$','interpreter','latex')
xlabel('$x$','interpreter','latex')
ylabel('$t$','interpreter','latex')
set(gca,'fontsize',16)

%f2
subplot(2,2,2)
surfl(real(f2));
shading interp
colormap (gray)
view(-20,60)
set(gca , 'YTick', numel(t)/4 * (0:4)),
set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca , 'XTick', linspace(1,numel(xi),3)),
set(gca , 'Xticklabel',{'-10', '0', '10'});
title('$f_2(x,t)$','interpreter','latex')
xlabel('$x$','interpreter','latex')
ylabel('$t$','interpreter','latex')
set(gca,'fontsize',16)

%sum of both function
subplot(2,2,3)
surfl(real(f));
shading interp
colormap (gray)
view(-20,60)
set(gca , 'YTick', numel(t)/4 * (0:4)),
set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca , 'XTick', linspace(1,numel(xi),3)),
set(gca , 'Xticklabel',{'-10', '0', '10'});
title('$f(x,t)=\bf{X}(x,t)$','interpreter','latex')
xlabel('$x$','interpreter','latex')
ylabel('$t$','interpreter','latex')
set(gca,'fontsize',16)

%the question is: if you measure (f), could you back up the two modes that
%build it, i.e., f1 and f2.

%% Perform DMD on data
%create DMD matrices
X1 = X(:,1:end-1);
X2 = X(:,2:end);

%STEP 1: SVD and rank-2 truncation
r = 2; %rank trauncation (2 modes)
[U,S,V] = svd(X1,'econ');
Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);

%Build Atilde and DMD modes: least square matrix/linear operator
Atilde = Ur'*X2*Vr/Sr;
%eigen decomposition
[W,D] = eig(Atilde); %eigenvectors and eigenvalues
%DMD modes 
Phi = X2*Vr/Sr*W; 

%DMD eigenvalues spectra
lambda = diag(D);
omega = log(lambda)/dt;

%reconstruction of modes/function in time
%compute DMD solution
x1 = X(:,1); %t=0
b = Phi\x1;
time_dynamics = zeros(r,length(t));
for iter = 1:length(t)
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
X_DMD = Phi * time_dynamics;

subplot(2,2,4)
surfl(real(X_DMD).')
shading interp
colormap gray
view (-20,60);
set(gca , 'YTick', numel(t)/4 * (0:4)),
set(gca , 'Yticklabel',{'0','\pi','2\pi','3\pi','4\pi'});
set(gca , 'XTick', linspace(1,numel(xi),3)),
set(gca , 'Xticklabel',{'-10', '0', '10'});
title('$\bf{X}_{\rm{DMD}}(x,t)$','interpreter','latex')
xlabel('$x$','interpreter','latex')
ylabel('$t$','interpreter','latex')
set(gca,'fontsize',16)

%check the eignevalues that all of them imaginary
% figure;
% plot(real(omega),imag(omega),'r*','Linewidth',2)
% xlabel('Real')
% ylabel('Imaginary')
% xlim([0 0.2])
% ylim([2 3])

%%Comparing PCA of the original dataset
%PCA
[u,s,v] = svd(X);
pc1 = u(:,1); %first PCA mode
pc2 = u(:,2); %second PCA mode
time_pc1 = v(:,1); %temporal evolution of pc1
time_pc2 = v(:,2); %temporal evolution of pc2

%visualise the comparison between PCA and DMD
figure;
%mode 1: spatial
subplot(1,2,1)
plot(xi,real(Phi(:,1)),'b','Linewidth',2)%DMD
hold on
plot(xi,real(pc1),'r','Linewidth',2)%PCA 
title('Mode 1')
axis square
xlabel('$x$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',16)

% %mode 2: spatial
subplot(1,2,2)
plot(xi,real(Phi(:,2)),'b','Linewidth',2)%DMD 
hold on
plot(xi,real(pc2),'r','Linewidth',2)%PCA
title('Mode 2')
axis square
xlabel('$x$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',16)

legend('DMD','PCA')
set(gca,'fontsize',16)

% %% svd, low rank structure
% [u,s,v] = svd(X);
% %eigenvalues of svd
% %plotting eigenvalues
% figure(2)
% plot(diag(s)/(sum(diag(s))),'ro') %it is obvious that we have two dominant modes
% %plotting eigenvectors
% figure(3)
% %plot the spatial modes
% subplot(2,1,1)
% plot(real(u(:,1:2)),'Linewidth',2)%modes
% %plot the temporal modes
% subplot(2,1,2)
% plot(real(v(:,1:2))) 





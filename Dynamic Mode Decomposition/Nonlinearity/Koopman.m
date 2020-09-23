clear all;
close all;

%%%%%%%%%%% NLS SOLVE - CREATION OF DATA
% space
L=30;
n=512;
x2=linspace(-L/2,L/2,n+1);
x=x2(1:n);
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
% time
slices=20;
t=linspace(0,2*pi,slices+1);
dt=t(2)-t(1);

%   nf=0.5;
%   slicesf=2*slices*nf;
%   tf=linspace(0,2*pi*nf,slicesf+1);

% initial conditions
N=2;
u=N*(sech(x)).';
ut=fft(u);
[t,utsol]=ode45('NLSE_solution_rhs',t,ut,[],k);
for j=1:length(t)
    usol(j,:)=ifft(utsol(j,:));  % bring back to space
end

subplot(1,5,1), waterfall(x,t,abs(usol)), colormap([0 0 0])
set(gca,'Ylim',[0 2*pi],'Ytick',[0 3 6],'Fontsize',14,'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('\xi','Fontsize',14), ylabel('time','Fontsize',14)
axis square

X = usol.';  % here is the data

%%%%%% body of DMD %%%%%%%%%%
X1 = X(:,1:end-1);
X2 = X(:,2:end);

[U2,Sigma2,V2] = svd(X1, 'econ');
r=10;
U=U2(:,1:r);
Sigma=Sigma2(1:r,1:r);
V=V2(:,1:r);
Atilde = U'*X2*V/Sigma;
[W,D] = eig(Atilde);
Phi=X2*V/Sigma*W;

mu=diag(D);
omega=log(mu)/dt;

y0 = Phi\u;  % pseudo-inverse initial conditions

u_modes = zeros(r,length(t));  % DMD reconstruction for every time point
for iter = 1:length(t)
    u_modes(:,iter) =(y0.*exp(omega*(t(iter))));
end
u_dmd = Phi*u_modes;   % DMD resconstruction with all modes

% figure(3), subplot(3,3,1), waterfall(x,1:r,abs(Phi).'), colormap([0 0 0])
% set(gca,'Xlim',[-8 8],'Xtick',[-8 0 8],'Ylim',[1 r],'Ytick',[1 r],'Fontsize',[14],'Zlim',[0 0.4],'Ztick',[0 0.2 0.4])
% figure(7), subplot(3,3,4), plot(omega,'ko','Linewidth',[2]), grid on, axis([-5 1 -20 20]), set(gca,'Fontsize',[14])
%
%
figure(1), subplot(1,5,2), waterfall(x,t,abs(u_dmd).'), colormap([0 0 0])
set(gca,'Ylim',[0 2*pi],'Ytick',[0 3 6],'Fontsize',14,'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('\xi','Fontsize',14), ylabel('time','Fontsize',14)
axis square

% figure(7), subplot(3,3,1), waterfall(x,tf,abs(u_dmd).'), colormap([0 0 0])
%    set(gca,'Ylim',[0 2*pi*nf],'Ytick',[0 3 6],'Fontsize',[14],'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
%
for j=1:length(t)
    error1(j)=norm(u_dmd(:,j)-usol(j,:).');
end


% Koopman

for jkoop=1:2
    if jkoop==1
        Y1=[X1; (X1.*abs(X1).^2)];
        Y2=[X2; (X2.*abs(X2).^2)];
    else if jkoop==2
            Y1=[X1; (abs(X1).^2)];
            Y2=[X2; (abs(X2).^2)];
        end
    end
    
    [U2,Sigma2,V2] = svd(Y1, 'econ');
    r=10;
    U=U2(:,1:r);
    Sigma=Sigma2(1:r,1:r);
    V=V2(:,1:r);
    Atilde = U'*Y2*V/Sigma;
    [W,D] = eig(Atilde);
    Phi2=Y2*V/Sigma*W;
    mu=diag(D);
    omega=log(mu)/dt;
    
    u2=[u; (u.*abs(u).^2)];
    y0 = Phi2\u2;  % pseudo-inverse initial conditions
    
    u_modes = zeros(r,length(t));  % DMD reconstruction for every time point
    for iter = 1:length(t)
        u_modes(:,iter) =(y0.*exp(omega*(t(iter))));
    end
    u_dmd2 = Phi2*u_modes;
    
    % DMD resconstruction with all modes
    u_dmd = u_dmd2(1:n,:);
    
    
    Phi = Phi2(1:n,:);
    
    if jkoop==1
        figure(1), subplot(1,5,2+jkoop), waterfall(x,t,abs(u_dmd).'), colormap([0 0 0])
        set(gca,'Ylim',[0 2*pi],'Ytick',[0 3 6],'Fontsize',14,'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
        axis square
    
    else if jkoop==2
        figure(1), subplot(1,5,2+jkoop), waterfall(x,t,abs(u_dmd).'), colormap([0 0 0])
        set(gca,'Ylim',[0 2*pi],'Ytick',[0 3 6],'Fontsize',14,'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
        axis square
        
        end
    end
end
    
    
    


clear all;
close all;

             %%% Define Discretised spatial and time variables %%% 

%Space
%-----
%Spatial variable on x direction
Nx=64; %number of spatial points
x2 = linspace(-10,10,Nx+1); %x variable
x = x2(1:Nx);
%Spatial variable on y direction
Ny=64; %number of spatial points
y2 = linspace(-10,10,Ny+1); %y variable
y = y2(1:Ny);
delta=20/Nx; %spatial step size
%Total matrix size
N = Nx*Ny;
%Time 
%-----
dt=0.3;  %time step size
tspan=0:dt:80; %time variable


                %%% Matrices of Numerical Derivatives %%%

e0=zeros(N,1); % vector of zeros
e1=ones(N,1); % vector of ones
e2=e1; % copy the one vector
e4=e0; % copy the zero vector
for j=1:Nx
    e2(Nx*j)=0; % overwrite every m^th value with zero
    e4(Nx*j)=1; % overwirte every m^th value with one
end
e3(2:N,1)=e2(1:N-1,1); 
e3(1,1)=e2(N,1); % shift to correct
e5(2:N,1)=e4(1:N-1,1); 
e5(1,1)=e4(N,1); % positions

%Create matrix A (Laplacian)
%---------------------------
% place diagonal elements
A = spdiags([e1 e1 e5 e2 -4*e1 e3 e4 e1 e1], ...
[-(N-Nx) -Nx -Nx+1 -1 0 1 Nx-1 Nx (N-Nx)],N,N);
A(1,1) = 2;
A = A/(delta^2);
%spy(A) %show the sparsity patterns of the matrix

%Create matrix B (first order derivative in x)
%---------------------------------------------
adds = (1/2/delta)*ones(N, 1);
subs = -1*adds;
B = spdiags([adds adds subs subs], [Nx -(N-Nx) (N-Nx) -Nx], N, N);
%spy(B) %show the sparsity patterns of the matrix

%Create matrix C (first order derivative in y)
%---------------------------------------------
C = spdiags([e5 -e2 e3 -e4], [-(Nx-1) -1 1 (Nx-1)], N, N);
C = (1/2/delta).*C;
%spy(C) %show the sparsity patterns of the matrix

                	     %%% Initial state %%%

%Create a meshgrid to let Matlab know the x and y directions
[X,Y] = meshgrid(x,y);
% Defining vorticity initial state/condition:
w0 = exp(-(X).^2 - ((Y.^2)/20)); %2D Gaussian Function 
w_vec = reshape(w0, [N 1]); %convert to a vector so it can be passed to ode45


                %%%Iterate and integrate over time %%%
%ode45
[Time,Omega] = ode45('fdm_rhs',tspan,w_vec,[],A,B,C); 


                       %%% plotting %%%
    
figure;
%t=0
subplot(3,3,1)
imagesc(abs(reshape(Omega(1,:),Nx,Nx)))
axis square
colormap jet
%t=0
subplot(3,3,2)
imagesc(abs(reshape(Omega(10,:),Nx,Nx)))
axis square
colormap jet
%t=0
subplot(3,3,3)
imagesc(abs(reshape(Omega(25,:),Nx,Nx)))
axis square
colormap jet
subplot(3,3,4)
%t=0
imagesc(abs(reshape(Omega(75,:),Nx,Nx)))
axis square
colormap jet
subplot(3,3,5)
%t=0
imagesc(abs(reshape(Omega(100,:),Nx,Nx)))
axis square
colormap jet
subplot(3,3,6)
%t=0
imagesc(abs(reshape(Omega(125,:),Nx,Nx)))
axis square
colormap jet
subplot(3,3,7)
%t=0
imagesc(abs(reshape(Omega(175,:),Nx,Nx)))
axis square
colormap jet
subplot(3,3,8)
%t=0
imagesc(abs(reshape(Omega(225,:),Nx,Nx)))
axis square
colormap jet
subplot(3,3,9)
%t=0
imagesc(abs(reshape(Omega(267,:),Nx,Nx)))
axis square
colormap jet


                       %%% Make a video %%%
%figure;
%loops = size(Omega,1);
%F(loops) = struct('cdata',[],'colormap',[]);
%for i = 1:loops
 %   test = Omega(i,:);
  %  test = reshape(test,64,64);
   % pcolor(abs(test));
    %axis off
    %colormap jet
    %shading interp
    %drawnow;
    %F(i) = getframe;
%end
%v = VideoWriter('Vorticity_FDM.avi');
%open(v);
%writeVideo(v, F);
%close(v);




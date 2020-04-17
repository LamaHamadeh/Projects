

close all;
clear all;

% Define Fourier points
N=31;
L=10;

% compute the Chebychev matrix
[D_unscaled, x_unscaled] = Cheb_Diff_Matrix(N-1);

%rescaling
D = (2/L).*D_unscaled;
x = (L/2).*x_unscaled;

%Imposig boundary conditions
D(N,:) = zeros(1,N);
D(1,:)   = zeros(1,N);

D(:,N) = zeros(N,1);
D(:,1)   = zeros(N,1);

% D = D(2:N,2:N);
% %or
% D = D(2:end-1,2:end-1);
% 
%second derivative
D2 = D^2; 
%first derivative
D1 = D;

%Define the second dimension
y=x;

% Build the Laplacian 2D matrix that satsfies boundary conditions
I = eye(length(D2));
A = kron(I,D2)+kron(D2,I);

%to calculate the d/dx
%the first derivative of x in 2D
B = kron(D1,I);

%to calculate the d/dy
%the first derivative of y in 2D
C = kron(I,D1);

%create the XY domain meshgrid
[X,Y]=meshgrid(x,y);

%initial condition
%vorticity
%option 1: %one Gaussian at the centre
w = exp(-2*X.^2-Y.^2/5); 

%option 2: %two same gaussians voticies next to each other which can be
%made to collide
%w = exp(-2*(X+2).^2-(Y).^2/20) + exp(-2*(X-2).^2-(Y).^2/20); 

%option 3: %multiple gaussians
% w = exp(-2*(X+6).^2-(Y).^2/10) + exp(-2*(X-6).^2-(Y).^2/10)...
%     +exp(-2*(X).^2-(Y+6).^2/10) + exp(-2*(X).^2-(Y-6).^2/10); 

%option 4: one positive Gaussian and one negative Gaussian
%w = exp(-2*(X+2).^2-(Y+6).^2/20) - exp(-2*(X-2).^2-(Y+6).^2/20); 

% plot the vorticity
% figure(1)
% pcolor(X,Y,real(w));
% shading interp
% colorbar
% colormap jet
% drawnow %to show the plot right away without waitigng for the code to be 
% finished running

%convert it ot a column
w = reshape(w,N^2,1); 

 
%solving diffusion equation
%time varibale
tspan=[0 4];

[Time,Omega] = ode45('cheb_rhs',tspan,w,[],A,B,C); %solve PDE

for j = 1:length(Time)
    pcolor(real(reshape(Omega(j,:),N,N)))
    shading interp
    drawnow
end

%                        %%% plotting %%%

% for j = 1:length(tspan)
%     w=ifft2(reshape(wt2sol(j,:),N+1,N+1));
%     figure(2)
%     subplot(4,4,j)
%     pcolor(x,y,real(w));
%     shading interp
%     colorbar
%     colormap jet
%     axis square
%     xlabel('x')
%     ylabel('y')
% end
% 

                       %%% Make a video %%%
% figure;
% loops = size(Sol,1);
% F(loops) = struct('cdata',[],'colormap',[]);
% for i = 1:loops
%     test = Sol(i,:);
%     test = reshape(test,N+1,N+1);
%     pcolor(abs(test));
%     %surfl(x,y,abs(test))
%     axis off
%     zlim ([0, 1]); %This is very important so we can lock the zoom-in behaviour
%     colormap jet
%     shading interp
%     drawnow;
%     F(i) = getframe;
% end
% v = VideoWriter('Heat_Cheb.avi');
% open(v);
% writeVideo(v, F);
% close(v);
% 
% 





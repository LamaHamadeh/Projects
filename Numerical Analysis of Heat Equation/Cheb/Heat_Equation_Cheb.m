

close all;
clear all;

% Define Fourier points
N=30;
% compute the Chebychev matrix
[D,x] = Cheb_Diff_Matrix(N);

%Boundary conditions
D(N+1,:) = zeros(1,N+1);
D(1,:)   = zeros(1,N+1);

%second derivative
D2 = D^2; 

%Define the second dimension
y=x;

% Build the Laplacian 2D matrix that satsfies boundary conditions
I = eye(length(D2));
L = kron(I,D2)+kron(D2,I);

%to calculate the d/dx+d/y 
%the first derivative in 2D, we use
%L = kron(I,I)+kron(D,I);


%create the XY domain meshgrid
[X,Y]=meshgrid(x,y);

%initial condition
U = exp(-(X.^2+Y.^2)/0.1);

%plotting
figure(1)
surfl(x,y,U)
shading interp
colormap jet
drawnow

%solving diffusion equation
u = reshape(U,(N+1)^2,1); %convert it ot a column
%time varibale
tspan=[0 0.1];

[Time,Sol] = ode23('cheb_rhs',tspan,u,[],L); %solve PDE


                       %%% Make a video %%%
figure;
loops = size(Sol,1);
F(loops) = struct('cdata',[],'colormap',[]);
for i = 1:loops
    test = Sol(i,:);
    test = reshape(test,N+1,N+1);
    pcolor(abs(test));
    %surfl(x,y,abs(test))
    axis off
    zlim ([0, 1]); %This is very important so we can lock the zoom-in behaviour on the z axis
    colormap jet
    shading interp
    drawnow;
    F(i) = getframe;
end
v = VideoWriter('Heat_Cheb.avi');
open(v);
writeVideo(v, F);
close(v);








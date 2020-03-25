
%Filtering is a general term for extracting information from a noisy signal. 
%Smoothing is a particular kind of filtering in which low-frequency 
%components are passed and high-frequency components are attenuated 
%(?low-pass filter?).


close all;
clear all;

%read image
A = imread('photo','jpeg'); %Glaxo Smith Kline Carbon Neutral Laboratory for Sustainable Chemistry
%resize image
%A = imresize(A,[350 350]); %if needed
%plotting
%subplot(1,2,1)
%imshow(A)
%axis square
%--------------------

%convert image to black and white
Abw = rgb2gray(A);
%convert the BW image from image format to a matrix
Abw= double(Abw);
%define number of pixels on each axis
[nx,ny]=size(Abw);
%add noise to the BW image
noise = randn(size(Abw,1),size(Abw,1));
B = Abw+20*noise;
%plotting
%subplot(1,2,2)
%imshow(uint8(B))
%colormap(gray(256))
%axis square
%--------------------

%Diffusion of an image/smooth the BW image

%Defining Numerical parameters:
%Spatial variable and its derivatives on x direction
x = linspace(0,1,nx); %spatial vector
dx = x(2) - x(1); %spatial step size
Ix = eye(nx); %unitary matrix
onex = ones(nx,1);
Dx = spdiags([onex -2*onex onex],[-1 0 1],nx,nx); %x 1D matrix

%Spatial variable on y direction
y = linspace(0,1,ny); %spatial vector
dy = y(2) - y(1); %spatial step size
Iy = eye(ny); %unitary matrix
oney = ones(ny,1);
Dy = spdiags([oney -2*oney oney],[-1 0 1],ny,ny); %y 1D matrix 

%2D Laplacian
L = kron(Dx,Iy)+kron(Ix,Dy); 

%initial conditions
U = Abw+20*noise; %our noisy BW image
U = reshape(U,nx*ny,1); %reshape as a vector

%Diffusion coeffieicnet
D = 0.005;
%Time variable
tspan = [0 0.002 0.004 0.006];

%solve diffusion equation
[Time,Sol] = ode45('Smoothing_rhs',tspan,U,[],L,D);

%--------------------

%plotting
%show the performance of the diffusion at certain time steps
for j = 1:length(tspan)
    Abw_smooth = uint8(reshape(Sol(j,:),nx,ny));
    subplot(2,2,j)
    imshow(Abw_smooth)
end
%--------------------


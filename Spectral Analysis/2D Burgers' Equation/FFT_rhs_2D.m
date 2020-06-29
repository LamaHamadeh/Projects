
function rhs = FFT_rhs_2D(tspan,U,dummy,Nx,Ny,Lx,Ly)

%Total matrix size
N = Nx * Ny;

%Fast Fourier Transform to the initial condition                 
Ut = fftshift(fft2(U)); %%FT+shift of the initial condition
Ut = reshape(Ut,N,1); %stack the data to a coloumn 
%-------------

%Define the 2D wave number
%x direction
kx = (2*pi/Lx)*[0:(Nx/2-1) (-Nx/2):-1]; 
kx(1) = 10^(-6);
kx = fftshift(kx);
%y diresction
ky = (2*pi/Ly)*[0:(Ny/2-1) (-Ny/2):-1]; 
ky(1) = 10^(-6);
ky = fftshift(ky);
%Create wave number mesh
[KX,KY] = meshgrid(kx,ky); %Nx X Ny grid
%convert to columns so they can pass to ode45
KX = reshape(KX,N,1);
KY = reshape(KY,N,1);
%--------------

%Advection Terms
%xaxis
duhatx = 1i .*KX .*Ut; %apply FT on the 1st derivative
dux = real(ifft2(ifftshift(duhatx))); %%inverse of FT
%yaxis
duhaty = 1i .*KY .*Ut; %apply FT on the 1st derivative
duy = real(ifft2(ifftshift(duhaty)));%inverse of FT

%Diffusion Terms
%xaxis
dduhatx = -(KX.^2) .*Ut; %apply FT on the 2nd derivative
ddux = real(ifft2(ifftshift(dduhatx))); %inverse of FT
%yaxis
dduhaty = -(KY.^2) .*Ut; %apply FT on the 2nd derivative
dduy = real(ifft2(ifftshift(dduhaty))); %inverse of FT
%-------------

%define Reynolds Number
Re  = 400; %not over 700 otherwise an error occurs: "Unable
%to meet integration tolerances"
nu = 1/Re;
%solve the right hand side
rhs = - U .* (dux + duy) + nu .* (ddux + dduy);
%the smaller the Reynolds Number, the smaller the shock wave produced.
%-------------

end

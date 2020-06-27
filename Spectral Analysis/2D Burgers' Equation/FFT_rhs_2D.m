
function rhs = FFT_rhs_1D(tspan,U,dummy,Nx,Ny,Lx,Ly)

%Total matrix size
N = (Nx * Ny);

%%% Fast Fourier Transform to the initial condition %%%                   
%FT+shift of the initial condition
Ut = fftshift(fft2(U));   
%stacked data to a coloumn to be passed later to ode45
Ut = reshape(Ut,N,1);

%x direction
kx = (2*pi/Lx)*[0:(Nx/2-1) (-Nx/2):-1]; 
kx(1) = 10^(-6);
%y diresction
ky = (2*pi/Ly)*[0:(Ny/2-1) (-Ny/2):-1]; 
ky(1) = 10^(-6);
%to give kx and ky the sense of direction
[KX,KY] = meshgrid(kx,ky); %N X N grid
%convert to columns so they can pass to ode45
KX = reshape(KX,N,1);
KY = reshape(KY,N,1);

%first derivative on x axis (advection)
duhatx = 1i .*KX .*Ut;
%inverse of FT
dux = ifft2(ifftshift(duhatx));
%first derivative on y axis (advection)
duhaty = 1i .*KY .*Ut;
%inverse of FT
duy = ifft2(ifftshift(duhaty));
%total advection on both axes
Advection = dux + duy;

%second derivative on the x axis (diffusion)
dduhatx = -(KX.^2) .*Ut;
%inverse of FT
ddux = ifft2(ifftshift(dduhatx));
%second derivative on the y axis (diffusion)
dduhaty = -(KY.^2) .*Ut;
%inverse of FT
dduy = ifft2(ifftshift(dduhaty));
%total diffusion on both axes
Diffusion = ddux + dduy;

%define Reynolds Number
Re  = 10;
%solve the right hand side
rhs = -U .*(Advection) + (1/Re) .*(Diffusion);
%the smaller the Reynolds Number, the smaller the shock wave produced.

end




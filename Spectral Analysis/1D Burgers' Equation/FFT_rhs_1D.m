
function rhs = FFT_rhs_1D(tspan,U,dummy,N,L)

%%% Fast Fourier Transform to the initial condition %%%                   
%FT+shift of the initial condition
Ut = fftshift(fft(U));   

%%% 1D Wave vector disretisation %%%                  
k = (2*pi/L)*[0:(N/2-1) (-N/2):-1]; %already shifted frequency domain
k(1) = 10^(-6);
k = fftshift(k);
k = reshape(k,N,1);

%first derivative (advection)
duhat = 1i .*k .*Ut;
%inverse of FT
du = ifft(ifftshift(duhat));

%second derivative (diffusion)
dduhat = -(k.^2) .*Ut;
%inverse of FT
ddu = ifft(ifftshift(dduhat));

%define diffusion coeffieicnt
%the larger the diffusion coeffiecient, the smaller the shock wave produced.
%i.e., less dynamics as the diffusion coefficient supresses the production
%of the shock wave.
diffusion  = 0.005;

%solve the right hand side
rhs = -U .*du + diffusion .*ddu;

end




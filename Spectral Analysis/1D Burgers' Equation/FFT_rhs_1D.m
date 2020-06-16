
function rhs = FFT_rhs_1D(tspan,U,dummy,N,L)

%%% Fast Fourier Transform to the initial condition %%%                   
%FT+shift of the initial condition
Ut = fft(U);   

%%% 1D Wave vector disretisation %%%                  
%option1:
k = (2*pi/L)*[0:(N/2-1) (-N/2):-1]; %already shifted frequency domain
k(1) = 10^(-6);
k = reshape(k,N,1);
%option2:
% k = (2*pi/L)*[-N/2:N/2-1]; %needs shifting
% k = fftshift(k);
% k = reshape(k,N,1);

%first derivative (advection)
duhat = 1i .*k .*Ut;
%inverse of FT
du = ifft(duhat);

%second derivative (diffusion)
dduhat = -(k.^2) .*Ut;
%inverse of FT
ddu = ifft(dduhat);

%define diffusion coeffieicnt
diffusion  = 0.01;
%solve the right hand side
rhs = -U .*du + diffusion .*ddu;

end




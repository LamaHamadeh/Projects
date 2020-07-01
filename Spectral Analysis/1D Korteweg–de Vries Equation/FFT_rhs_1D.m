
function rhs = FFT_rhs_1D(tspan,U,dummy,N,L)


%Fast Fourier Transform to the initial condition                 
Ut = fftshift(fft(U));   
% Ut = reshape(Ut,N,1);

%1D Wave vector disretisation              
k = (2*pi/L)*[0:(N/2-1) (-N/2):-1];
k(1) = 10^(-6);
% k = (2*pi/L)*(-N/2:N/2-1)';
k = fftshift(k);
k = reshape(k,N,1);

%first derivative (advection)
duhat = 1i *k .*Ut;
du = real(ifft(fftshift(duhat))); %inverse of FT

%third derivative (diffusion)
ddduhat = -1i * (k.^3) .*Ut;
dddu = real(ifft(fftshift(ddduhat))); %inverse of FT

%define diffusion coeffieicnt
diffusion  = 1;
advection = 6;

%solve the right hand side
rhs = - advection .* U .* du + diffusion .* dddu;

end



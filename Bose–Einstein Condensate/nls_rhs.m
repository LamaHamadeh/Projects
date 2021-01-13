

%i*ut+(1/2)*uxx+u.*|u|^2=0

function rhs = nls_rhs(t, ut, dummy, k)
%inverse of Fourier Transform of u
u = ifft(ut);
%solve the right hand side
rhs = -(1i/2)*(k.^2).*ut+1i*fft((abs(u).^2).*u);


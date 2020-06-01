clear all;
close all;

N = 101;
L = 8*pi;
dx = L/(N-1);
x = -L/2:dx:L/2-dx;
f = cos(x).*exp(-x.^2/25);
df = -(sin(x).*exp(-x.^2/25) + (2/25)*x.*cos(x).*exp(-x.^2/25));
Nx = max(size(f));
plot(x,f,'b',x,df,'k','LineWidth',1.5)
hold on

% Approximate derivative using finite Difference...
for k=1:length(f)-1
    dfFD(k) = (f(k+1)-f(k))/dx; %Forward difference
end
dfFD(end+1) = dfFD(end);
plot(x,dfFD,'r--')

% Derivative using FFT (spectral derivative)
fhat = fft(f);
k = (2*pi/L)*[-Nx/2:Nx/2-1];
k = fftshift(k);  % important because fft has weird ordering
dfhat = 1i*k.*fhat;
dfspec = real(ifft(dfhat));
plot(x,dfspec,'go')
xlabel('Spatial variable, x')
ylabel('Derivative of Function')

% Plotting commands
axis([-L/2 L/2 -1 1])
legend('Function','True Derivative','Finite Difference','FFT Derivative')

%reason behnind FD is less accurat than FFT
% small neighbour of points around where I'm taking the derivative
% and I am estimating what the curvature (the derivative) is around
% these piints. 
% If I use larger and larger radius of points, I'll get more and more
% accurate shceme. 
%However, in the spectral scheme, FFT, for each frequncy I have a 
% sine wave that covers all of the points in my domain. So doing 
% the derivative in the frequency domain is like takign a stencil that
%takes all the points in my domain, i.e., infifinte stencil.



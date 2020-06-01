clear all;
close all;

%create a spatial grid
N = 1001;
L = 100;
dx = L/(N);
x = -L/2:dx:L/2;
f = 0*x; %initial temperature distribution
f((L/2 - L/10)/dx:(L/2 + L/10)/dx) = 1; %a square wave
plot(x,f)

a = 1; %thermal diffusion constant
y = x; %summy variable for the convolution integral

dt = 0.1;
fhat = fft(f);

for k = 1:100
    t = k*dt; 
    Nx = max(size(f));
    kappa = (2*pi/L)*[-Nx/2:Nx/2-1];
    kappa = fftshift(kappa);  % important because fft has weird ordering
    uhat = fhat.*exp(-(a^2)*t*kappa.^2);
    u = ifft(uhat);
    plot(x,u,'k')
    xlabel('Spatial variable, x')
    ylabel('Temperature, u(x,t)')
    title(['Time, t=',num2str(t)])
    axis([-10 10 -.1 1.1])

    pause(0.1)
end
clear all;
close all;

%create a spatial grid
N = 1001;
L = 100;
dx = L/(N);
x = -L/2:dx:L/2;
f = 0*x; %initial temperature distribution
f((L/2 - L/10)/dx:(L/2 + L/10)/dx) = 1; %a square wave
fhat = fft(f);
a = 5;
y = x;
dt = 0.1;
for k = 1:100
    t = k*dt;
    for j = 1:length(x)
        xi = x(j) - y;
        gxy = (1/(2*a*sqrt(pi*t)))*exp(-xi.^2/(4*a*a*t));
        u(j) = dot(gxy,f)*dx;
    end
    plot(x,u,'k')
    hold on
    Nx = max(size(f));
    kappa = (2*pi/L)*[-Nx/2:Nx/2-1];
    kappa = fftshift(kappa);  % important because fft has weird ordering
    uhat = fhat.*exp(-(a^2)*t*kappa.^2);
    u = ifft(uhat);
    plot(x,u,'r--')
    xlabel('Spatial variable, x')
    ylabel('Temperature, u(x,t)')
    title(['Time, t=',num2str(t)])
    hold off
    axis([-L/2 L/2 -.1 1.1])
    legend('Analytical','Numerical') 
    % to make them agree, we need to make the domain even longer
    pause(0.1)
end


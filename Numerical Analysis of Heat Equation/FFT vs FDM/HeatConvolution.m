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

a = 3; %thermal diffusion constant
y = x; %summy variable for the convolution integral

dt = 0.1;
for k = 1:100 %100 time steps
    t = k*dt;
    for j = 1:length(x) %compute the convolution integral at each x point
        xi = x(j) - y;
        gxy = (1/(2*a*sqrt(pi*t)))*exp(-xi.^2/(4*a*a*t));
        u(j) = dot(gxy,f)*dx; %temperature at the jth spatial locatiob
    end
    plot(x,u,'k')
    xlabel('Spatial variable, x')
    ylabel('Temperature, u(x,t)')
    title(['Time, t=',num2str(t)])    
    axis([-L/2 L/2 -.1 1.1])

    pause(0.1)
end

clear all;
close all;


%Spatial variable on x direction
L=100; %domain
xmin=-L/2; %minimum boundary
xmax=L/2; %maximum boundary 
dx=0.05; %spatial step size
N=(xmax-xmin)/dx; %number of spatial points
x=linspace(xmin,xmax,N); %spatial vector
%initialise temeprature distribution
f = 0*x; 
%initial condition: square function
f((L/2 - L/10)/dx:(L/2 + L/10)/dx) = 1; 
%diffusion coeffieicnt
a = 10; 
%define the diffusion time step
dt = 0.1;
%define a dummy variable for the convolution integral
y = x;

Time = 100;

for i = 1:Time
    
    t = i*dt; %increase by one time step
    
    %Analytical: convolution integral
    for j = 1:length(x) %compute the convolution integral at each x point
        xi = x(j) - y;
        gxy = (1/(2*a*sqrt(pi*t)))*exp(-xi.^2/(4*a*a*t));
        u(j) = dot(gxy,f)*dx; %temperature at the jth spatial location
    end
    %plotting
    plot(x,u,'k')
    hold on
    
    %Numerical: Fourier transform
    %method 1
    %Fast Fourier transform of the function f
    fhat = fft(f);
    %Define the wave vector
    k = (2*pi/L)*[-N/2:N/2-1]; 
    k = fftshift(k);
    uhat = fhat.*exp(-(a^2)*t*k.^2);
    %inverse of FT to the shifted function
    u = ifft(uhat);
    %plotting
    plot(x,u,'r--')

%     %method 2
%     %Fast Fourier transform of the function f
%     fhat = fft(f);
%      %Define the wave vector
%     k = (2*pi/L)*[0:(N/2-1) (-N/2):-1]; %different than the previously defined k
%     %k(1) = 10^(-6);
%     uhat = fftshift(fhat.*exp(-(a^2)*t*k.^2));   
%     %inverse of FT to the shifted function
%     u = ifft(uhat);
%     %plotting
%     plot(x,abs(u),'r--')

    xlabel('Spatial variable, x')
    ylabel('Temperature, u(x,t)')
    title(['Time, t=',num2str(t)])
    hold off
    axis([-L/2 L/2 -.1 1.1])
    legend('Analytical','Numerical')
    
    pause(0.1)
end

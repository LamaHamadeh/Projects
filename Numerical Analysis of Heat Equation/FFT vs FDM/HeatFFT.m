clear all;
close all;

%Spatial variable on x direction
L=20; %domain
xmin=-L/2; %minimum boundary
xmax=L/2; %maximum boundary 
dx=0.05; %spatial step size
N=(xmax-xmin)/dx; %number of spatial points
x=linspace(xmin,xmax,N); %spatial vector

%initialise temeprature distribution
f = 0*x; 

%initial condition: square function
f((L/2 - L/10)/dx:(L/2 + L/10)/dx) = 1; 

%plotting the initial condition
plot(x,f)

%diffusion coeffieicnt
a = 1; 

%define the diffusion time step
dt = 0.1;

%Define the wave vector
k = (2*pi/L)*[0:(N/2-1) (-N/2):-1]; 
%k(1) = 10^(-6);

%solve the 1D diffusion equation as a function of time.
Time = 100;

for i = 0:Time     
    t = i*dt; %increase by one time step   
    %Fast Fourier transform of the function f
    fhat = fft(f);
    uhat = fftshift(fhat.*exp(-(a^2)*t*k.^2));    
    %inverse of FT to the shifted function
    u = ifft(uhat);    
    %plotting
    plot(x,abs(u),'k') %the absolute value of the function
    xlabel('Spatial variable, x')
    ylabel('Temperature, u(x,t)')
    title(['Time, t=',num2str(t)])
    axis([-15 15 -.1 1.1])
    pause(0.1)
end


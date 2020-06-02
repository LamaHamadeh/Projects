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

%define a dummy variable for the convolution integral
y = x;

for k = 1:100 %100 time steps
    t = k*dt; %increase by one time step    
    for j = 1:length(x) %compute the convolution integral at each x point
        xi = x(j) - y;
        gxy = (1/(2*a*sqrt(pi*t)))*exp(-xi.^2/(4*a*a*t));
        u(j) = dot(gxy,f)*dx; %temperature at the jth spatial location
    end
    %plotting
    plot(x,u,'k')
    xlabel('Spatial variable, x')
    ylabel('Temperature, u(x,t)')
    title(['Time, t=',num2str(t)])    
    axis([-L/2 L/2 -.1 1.1])
    pause(0.1)
end

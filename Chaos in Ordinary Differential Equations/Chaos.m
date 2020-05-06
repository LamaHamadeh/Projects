

clear all; 
close all;

%define time domain
tmin = 0;
dt = 0.01;
tmax = 10;
tspan = tmin:dt:tmax;

%define initial conditions for x,y,z
x0 = [-8; 8; 27];

%define Lorentz parameters
sigma = 10; 
b=8/3; 
r = 28;

%solve the three ODEs using ode45
[t,Y1] = ode45('Lorentz', tspan, x0, [], sigma, b, r);
%define dimensions
x1 = Y1(:,1);
y1 = Y1(:,2);
z1 = Y1(:,3);
%plotting
%3D
figure(1)
plot3(x1,y1,z1,'b', 'LineWidth',2), grid on
xlabel('x')
ylabel('y')
zlabel('z')
legend('x0')

hold on

%define perturbation/butterfly effect
epsilon = 0.01*ones(3,1);
xperturbed = x0 + epsilon;

%solve the three ODEs using ode45
[t,Y2] = ode45('Lorentz', tspan, xperturbed, [], sigma, b, r);
%define dimensions
x2 = Y2(:,1);
y2 = Y2(:,2);
z2 = Y2(:,3);
%plotting
%3D
figure(1)
plot3(x2,y2,z2,'r', 'LineWidth',2), grid on
xlabel('x')
ylabel('y')
zlabel('z')
legend('x0','xperturbed')

%2D
figure(2)
subplot(3,1,1), plot(t,x1,'b',t,x2,'bo')
xlabel('t')
ylabel('x')
subplot(3,1,2), plot(t,y1,'r',t,y2,'ro')
xlabel('t')
ylabel('y')
subplot(3,1,3), plot(t,z1,'g',t,z2,'go')
xlabel('t')
ylabel('z')


%analyse the difference between x0 and xperturbed
for i = 1:length(tspan)
    diff(i) = norm(Y2(i,:) - Y1(i,:));
end
figure(3)
plot(diff)
xlabel('Time')
ylabel('Difference between x0 and xperturbed')

%We can see that at the very beginning of time, they show exponential 
%divergence of the two trajectories. this little epsilon is multiplied 
%by an exponential function at the early time. This is another hallmark
%for chaotic systems/functions.



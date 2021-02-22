
clear all;
close all;

%Define spring-mass-damper system parameters
w = 2*pi; %natural frequency
d = 0.25; %Damping ratio

%Spring Mass damper System matrix
A = [0 1 ; -w^2 -2*d*w];

%Define time variable
dt = 0.01; %time step
T = 10; %total time domain of integration

%initial condition
x0 = [2 ; 0]; %[x = 2, v = 0].

%Iterate Forward-Euler Scheme
xF(:,1) = x0; %state
tF(1) = 0; %time

for k = 1:T/dt
    tF(k+1) = k+dt;
    xF(:,k+1) = (eye(2) + dt*A)*xF(:,k);
end

%Iterate Backward-Euler Scheme
xB(:,1) = x0; %state
tB(1) = 0; %time

for k = 1:T/dt
    tB(k+1) = k+dt;
    xB(:,k+1) = inv(eye(2) - dt*A)*xB(:,k);
end

%4th-order Runge-Kutta algorithm/integrator
[t,xGood] = ode45(@(t,x) A*x, 0:dt:T, x0); %we are using function hamdles

%Visualisation
figure;
subplot(1,2,1)
time = linspace(0,T,1001); %construct a time variable for the x axis
plot(time,xF(1,:),'b')
hold on
plot(time,xB(1,:),'r')
hold on
plot(time,xGood(:,1),'k')
xlabel('Time')
ylabel('Position')
legend('Forward-Euler','Backward-Euler','Runge-Kutta')
axis square

%phase space diagram
subplot(1,2,2)
plot(xF(1,:),xF(2,:),'b')
hold on
plot(xB(1,:),xB(2,:),'r')
hold on
plot(xGood(:,1),xGood(:,2),'k')
xlabel('Position (x)')
ylabel('Velocity (v)')
legend('Forward-Euler','Backward-Euler','Runge-Kutta')
axis square



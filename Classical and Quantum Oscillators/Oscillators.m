

clear all; %clear all previously defined variables
close all; %close all previosly defined figures


% %%%Damped Harmonic Oscillator
% %----------------------------
% %In classical mechanics, a harmonic oscillator is a system that, 
% %when displaced from its equilibrium position, experiences a restoring 
% %force F proportional to the displacement x.
% 
% %If F is the only force acting on the system, the system is called a simple
% %harmonic oscillator, and it undergoes simple harmonic motion: sinusoidal
% %oscillations about the equilibrium point, with a constant amplitude and 
% %a constant frequency (which does not depend on the amplitude).
% 
% %If a frictional force (damping) proportional to the velocity is also
% %present, the harmonic oscillator is described as a damped oscillator.
% 
% %define natural frequency 
% w = 2*pi;
% 
% %define damping ratio
% d = 0.50;
% 
% %initial conditions of the position and its first drivative
% y0 = [2,0];
% 
% %Time variable
% tmin=0; %minimum time
% tmax=10; %maximum time
% tspan = [tmin tmax];
% 
% %solve ODE using explicit Runge-Kutta (4,5)
% [Time,Position] = ode45('Oscillators_rhs',tspan,y0,[],d,w);
% 
% %plotting
% plot(Time,Position(:,1),'-o',Time,Position(:,2),'-o')
% title('Solution of Damped Harmonic Oscillator with ODE45');
% xlabel('Time');
% ylabel('Solution');
% legend('Position','Velocity')


% %%%Van Der Pol Equation
% %----------------------
% %Definition:
% %In dynamics, the Van der Pol oscillator is a non-conservative oscillator 
% %with non-linear damping. It evolves in time according to a second-order 
% %ordinary differential equation
% 
% %define the  scalar parameter indicating the nonlinearity and the strength 
% %of the damping.
% mu = 1;
% 
% %initial conditions of the position and its first drivative
% y0 = [2,0];
% 
% %Time variable
% tmin=0; %minimum time
% tmax=20; %maximum time
% tspan = [tmin tmax];
% 
% %solve ODE using explicit Runge-Kutta (4,5)
% [Time,Position] = ode45('Oscillators_rhs',tspan,y0,[],mu);
% 
% %plotting
% plot(Time,Position(:,1),'-o',Time,Position(:,2),'-o')
% title('Solution of van der Pol Equation (\mu = 1) with ODE45');
% xlabel('Time');
% ylabel('Solution');
% legend('Position','Velocity')


%%%Quantum Harmonic Oscillator 
%-----------------------------
%Definition: 
%The quantum harmonic oscillator is the quantum-mechanical 
%analog of the classical harmonic oscillator. Because an arbitrary smooth 
%potential can usually be approximated as a harmonic potential at the 
%vicinity of a stable equilibrium point, it is one of the most important 
%model systems in quantum mechanics. Furthermore, it is one of the few 
%quantum-mechanical systems for which an exact, analytical solution is 
%known.

%Spatial variable
L=8; %spatial domain
dq=0.05; %spatial step size
qmin=-L; %minimum boundary
qmax=L; %maximum boundary 
N=(qmax-qmin)/dq; %number of spatial points
q=linspace(qmin,qmax,N); %spatial vector

%Define the harmonic oscillator length scale
alpha = 0.255; 

% Constructing the spatial matrix using Finite-Difference Method (FDM)
A=zeros(N,N); %matrix of zero elements
%the diagonal elements
for n=1:N  %the number of colomns/rows
    A(n,n) = 2/dq^2 + alpha * q(n)^2;  %the value of each diagonal element 
end

%the off-diagonal elements
for n=1:N-1       
    A(n+1,n) = -1/dq^2; %the value of each lower off-diagonal element
end 
for n=2:N        
    A(n-1,n) = -1/dq^2; %the value of each upper off-diagonal element
end

%computing eigenfunctions and eigenvalues
[Eigen_Fun,Eigen_Val] = eig(A);
D=diag(Eigen_Val); %eigenvalues

%normalisation
%1st mode
norm1 = trapz(q,Eigen_Fun(:,1).*Eigen_Fun(:,1));
Eigen_Fun1 = Eigen_Fun(:,1)/sqrt(norm1);
%2nd mode
norm2 = trapz(q,Eigen_Fun(:,2).*Eigen_Fun(:,2));
Eigen_Fun2 = Eigen_Fun(:,2)/sqrt(norm2);
%3rd mode
norm3 = trapz(q,Eigen_Fun(:,3).*Eigen_Fun(:,3));
Eigen_Fun3 = Eigen_Fun(:,3)/sqrt(norm3);
%4th mode
norm4 = trapz(q,Eigen_Fun(:,4).*Eigen_Fun(:,4));
Eigen_Fun4 = Eigen_Fun(:,4)/sqrt(norm4);

%plotting
%First four normalised eigenfunctions 
figure(1)
plot(q,Eigen_Fun1,'b',q,Eigen_Fun2,'ro',...
    q,Eigen_Fun3,'g--',q,Eigen_Fun4,'m.','linewidth',2)
xlabel('x')
ylabel('Normlaised Eigenfuncions \psi(x)')
legend('1st','2nd','3rd','4th')
%the probability density distribution of the quantum state n=12
figure(2)
norm12 = trapz(q,Eigen_Fun(:,12).*Eigen_Fun(:,12));
Eigen_Fun12 = Eigen_Fun(:,12)/sqrt(norm12);
plot(q,abs(Eigen_Fun12).^2,'b')
xlabel('x')
ylabel({'The probability density distribution for finding ';...
    'the quantum harmonic oscillator in its n=12 quantum state: |\psi(x)_{12}|^2'})

%A function to caluclate the time-dependent temperature distribution (TDTD)

%% Define the function
%--------------------
function dT=dTDistribution(tspan,T0,dummy,A,dx) 
%the order after 'dummy'is important and should match with ode45 syntax

dT = A*T0; 
end




%A function to caluclate the time-dependent temperature distribution (TDTD)

%% Define the function
%--------------------
function dT=dTDistribution_3clouds(tspan,U,dummy,A,delta,D) 
%the order after 'dummy'is important and should match with ode45 syntax

dT = D*A*U; 
end

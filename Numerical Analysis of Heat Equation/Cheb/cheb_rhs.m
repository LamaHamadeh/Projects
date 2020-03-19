

function rhs=cheb_rhs(tspan,U,dummy,L) 
%the order after 'dummy'is important and should match with ode45 syntax

rhs = L * U; 
end
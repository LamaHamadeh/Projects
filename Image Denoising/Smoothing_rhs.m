
%define the function
function rhs = Smoothing_rhs(tspan,U,dummy,L,D)

%compute the right hand side of the diffusion equation
rhs = D*L*U;

end
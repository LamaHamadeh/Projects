

function rhs=cheb_rhs(tspan,w,dummy,A,B,C) 
%the order after 'dummy'is important and should match with ode45 syntax

%Define viscosity
visocity = 0.001;

%compute the streamfunction at t=0
psi0 = A\w;

%compute the matrices multiplications
psix   = B * psi0;
psiy   = C * psi0;
wx     = B * w;
wy     = C * w;
diff_w = A * w;

%compute vorticity at a future step and iterate/loop over all time steps
rhs = visocity.*diff_w - psix.*wy + psiy.*wx;

end
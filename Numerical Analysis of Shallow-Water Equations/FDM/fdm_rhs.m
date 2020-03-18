

function rhs = fdm_rhs(tspan,w_vec,dummy,A,B,C) 

%Define viscosity
visocity = 0.001;

%compute the streamfunction at t=0
psi0 = A\w_vec;

%compute the matrices multiplications
psix   = B * psi0;
psiy   = C * psi0;
wx     = B * w_vec;
wy     = C * w_vec;
diff_w = A * w_vec;

%compute vorticity at a future step and iterate/loop over all time steps
rhs = visocity.*diff_w - psix.*wy + psiy.*wx;

end
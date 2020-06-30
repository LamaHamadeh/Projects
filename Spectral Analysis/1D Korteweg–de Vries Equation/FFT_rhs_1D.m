
function rhs = FFT_rhs_1D(tspan,U,dummy,du,dddu)


%define diffusion coeffieicnt
diffusion  = 0.0005;
% diffusion = (0.022)^2;

%solve the right hand side
rhs = - U .* du + diffusion .* dddu;

end




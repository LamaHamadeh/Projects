
function rhs = FFT_rhs(tspan,Ut,dummy, k)

%solve the right hand side
rhs = 1i.*k.*Ut;

end




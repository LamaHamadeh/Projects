
function rhs = FFT_rhs(tspan,Ut,dummy, k)

%diffusion coeffieicnt
D = 0.01;
%solve the right hand side
rhs = -D.*k.^2.*Ut;

end




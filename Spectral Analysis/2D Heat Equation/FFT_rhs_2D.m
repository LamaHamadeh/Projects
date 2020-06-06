
function rhs = FFT_rhs_2D(tspan,Ut,dummy, KX, KY)

%diffusion coeffieicnt
D = 0.01;
%solve the right hand side
rhs = -D.*(KX.^2+KY.^2).*Ut;

end





function rhs = FFT_rhs(tspan,Ut2,dummy, KX, KY,D)

rhs = -D.*(KX.^2+KY.^2).*Ut2;

end




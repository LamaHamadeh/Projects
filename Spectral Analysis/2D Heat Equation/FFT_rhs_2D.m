
function rhs = FFT_rhs_2D(tspan,Ut,dummy, KX, KY)

%diffusion coeffieicnt
D = 0.01;

rhs = -D.*(KX.^2+KY.^2).*Ut;

end




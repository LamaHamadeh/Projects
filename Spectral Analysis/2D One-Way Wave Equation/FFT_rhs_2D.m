function rhs = FFT_rhs_2D(tspan,Ut,dummy, KX, KY)

%solve the right hand side
%rhs = -1i.*(KX+KY).*Ut;

%To make the wave go to the other way, the right hand side can be put as:
rhs = 1i.*(KX+KY).*Ut;

end


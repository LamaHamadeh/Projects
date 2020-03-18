
function rhs = spc_rhs(tspan,wt2,dummy,Kderv, KX, KY, nu,nx,ny,N)

%reshape the Fourier version of the intiial state
wt = reshape(wt2,nx,ny);

%calculate the streamfunction
psi = -wt./Kderv;

%compute the matrices multiplications
psix = real(ifft2(1i*KX.*psi));
psiy = real(ifft2(1i*KY.*psi));
wx = real(ifft2(1i*KX.*wt));
wy = real(ifft2(1i*KY.*wt));

%compute vorticity at a future step and iterate/loop over all time steps
rhs = -nu*Kderv.*wt+fft2(psiy.*wx-psix.*wy);

%reshape the solution as a vector
rhs = reshape(rhs,N,1);







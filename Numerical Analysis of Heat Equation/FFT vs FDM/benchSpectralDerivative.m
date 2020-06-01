clear all;
close all;

errFD = []; %FD eror
errFFT = []; %FFT error
Nvec = [8 16 32 64 128 256 512 1024 2048 4096]; %number of points
for N=Nvec
    L = 8*pi;
    dx = L/N;
    x = -L/2:dx:L/2-dx;
    f = cos(x).*exp(-x.^2/25);
    df = -(sin(x).*exp(-x.^2/25) + (2/25)*x.*cos(x).*exp(-x.^2/25));
    Nx = max(size(f));        
    plot(x,df,'k','LineWidth',1.5)
    hold on
    
    % Approximate derivative using finite Difference...
    for k=1:length(df)-1
        dfFD(k) = (f(k+1)-f(k))/dx;
    end
    dfFD(end+1) = dfFD(end);
    errFD = [errFD; dfFD((Nx)/2)-df((Nx)/2)];
    plot(x,dfFD,'b--')
    
    % Derivative using FFT (spectral derivative)
    fhat = fft(f);
    kappa = (2*pi/L)*[-Nx/2:Nx/2-1];
    kappa = fftshift(kappa);  % important because fft has weird ordering

    dfhat = 1i*kappa.*fhat;
    dfFFT = real(ifft(dfhat));
    errFFT = [errFFT; dfFFT((Nx)/2)-df((Nx)/2)];
    plot(x,dfFFT,'r--')
    
    % Plotting commands
    axis([-L/2 L/2 -1 1])
    legend('True Derivative','Finite Difference','FFT Derivative')
    xlabel('Spatial variable, x')
    ylabel('Derivative of Function')
    hold off
end

%
figure
loglog(Nvec,abs(errFD),'b-*')
hold on
loglog(Nvec,abs(errFFT),'r-*')
xlabel('Number of Discretization Points, N')
ylabel('Error')
legend('Finite Difference','Spectral Derivative')
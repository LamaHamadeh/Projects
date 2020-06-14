close all;
clear all;

%%% Wave number and wavelength
%Define the space domain
x = 0:0.01:6*pi;
%plotting different signals with different wave numbers
figure;
%k=1
y = cos(x);
subplot(2,2,1)
plot(x,y,'b','LineWidth',2)
xlim([0 6*pi])
xticks(0:pi:6*pi)
xticklabels({'0','\pi','2\pi','3\pi','4\pi','5\pi','6\pi'})
set(gca,'FontSize',16)
title('$k=1,\lambda = 2\pi$','Interpreter','latex')
%k=2
y = cos(2*x);
subplot(2,2,2)
plot(x,y,'b','LineWidth',2)
xlim([0 6*pi])
xticks(0:pi:6*pi)
xticklabels({'0','\pi','2\pi','3\pi','4\pi','5\pi','6\pi'})
set(gca,'FontSize',16)
title('$k=2,\lambda = \pi$','Interpreter','latex')
%k=3
y = cos(3*x);
subplot(2,2,3)
plot(x,y,'b','LineWidth',2)
xlim([0 6*pi])
xticks(0:pi:6*pi)
xticklabels({'0','\pi','2\pi','3\pi','4\pi','5\pi','6\pi'})
set(gca,'FontSize',16)
title('$k=3,\lambda = \frac{2\pi}{3}$','Interpreter','latex')
%k=4
y = cos(4*x);
subplot(2,2,4)
plot(x,y,'b','LineWidth',2)
xlim([0 6*pi])
xticks(0:pi:6*pi)
xticklabels({'0','\pi','2\pi','3\pi','4\pi','5\pi','6\pi'})
set(gca,'FontSize',16)
title('$k=4,\lambda = \frac{\pi}{2}$','Interpreter','latex')

h = suptitle('$f(x) = \cos(kx)$');
set(h,'FontSize',26,'FontWeight','normal','Interpreter','latex')

%--------------------------------------------------------------------------

%%% spectrogram of different wavenumber of cosine waves.
%Space analysis
%---------------
%Space domain
L = 6*pi; %is the range of the space domain
%Define space variable
n = 1884; %number of points in space domain
x2 = linspace(0,L,n+1);
x = x2(1:n); %space domain
%Define signal
wave_No = 7; %wave number
S = cos(wave_No*x); 
%plotting the signal
figure;
subplot(2,1,1)
plot(x,S,'b')
title('Signal in space domain')
xticks(0:pi:6*pi)
xticklabels({'0','\pi','2\pi','3\pi','4\pi','5\pi','6\pi'})
xlabel('Space')
ylabel('Signal')
set(gca,'FontSize',16)

hold on

%Frequency analysis
%-----------------
%Define wave number vector
k = (2*pi./L).*[0:n/2-1 -n/2:-1];
%Fourier shift to the wave number
ks = fftshift(k); 
%Fourier Transform of the signal
St = fft(S);
%plotting
subplot(2,1,2)
plot(ks,abs(fftshift(St))/max(abs(fftshift(St))),'r'); %normalised frequency intensity
title('Signal in frequency domain')
xlabel('Wave number ($m^{-1}$)','Interpreter','latex')
ylabel('Intensity')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
xlim([-10 10])

% Time frequency analysis (spectrogram)
% ----------------------------------------

%Initialise spectrogram matrix
Sgt_spec = []; 
%Define the time variable for iteration
xslide = 0:0.1:L; 
%Loop over all time steps
for j=1:length(xslide)
    width = 1;
    g = exp(-width*(x-xslide(j)).^2);
    
    %Gabor Signal
    Sg = g.*S; 

    %Fourier Transofrm of the resulting signal
    Sgt = fft(Sg);

    %store the spectrogram data at each iteration
    Sgt_spec = [Sgt_spec;abs(fftshift(Sgt))/max(abs(fftshift(Sgt)))]; 

    %plotting
    figure(2)
    %the signal and the sliding window function
    subplot(3,1,1)
    plot(x,S,'k',x,g,'r')
    title('Signal with Gabor function in space domain')
    xlabel('Space')
    xticks(0:pi:6*pi)
    xticklabels({'0','$\pi$','$2\pi$','$3\pi$','$4\pi$','$5\pi$','$6\pi$'})
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',16)
    %the convolution of the signal with the window function
    subplot(3,1,2)
    plot(x,Sg,'k')
    title('Convolution of the Signal with Gabor Transform in space domain')
    xticks(0:pi:6*pi)
    xticklabels({'0','$\pi$','$2\pi$','$3\pi$','$4\pi$','$5\pi$','$6\pi$'})
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',16)
    ylim([-1 1])
    %the normalised resulting signal in frequency domain
    subplot(3,1,3)
    plot(ks,abs(fftshift(Sgt))/max(abs(fftshift(Sgt))))
    xlabel('Wave number ($m^{-1}$)','Interpreter','latex')
    title('Fourier Transofrm of the new signal in frequency domain')
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',16)
    xlim([-10 10])
    drawnow

    pause(0.1)
end

%plotting the spectrogram
figure;
subplot(2,1,1)
plot(x,S,'b')
title('Signal in space domain')
xticks(0:pi:6*pi)
xticklabels({'0','\pi','2\pi','3\pi','4\pi','5\pi','6\pi'})
xlabel('Space')
ylabel('Signal')
xlim([0 6*pi])
set(gca,'FontSize',16)
subplot(2,1,2)
pcolor(xslide,ks,Sgt_spec.') 
shading interp
set(gca, 'Ylim', [-10 10])
colormap(jet)
xlabel('Space')
ylabel('Wave number ($m^{-1}$)','Interpreter','latex')
title('Amplitude of a particular wave number at a specific space location')
xticks(0:pi:6*pi)
xticklabels({'0','$\pi$','$2\pi$','$3\pi$','$4\pi$','$5\pi$','$6\pi$'})
set(gca,'TickLabelInterpreter','latex')
xlim([0 6*pi])
set(gca,'FontSize',16)

%--------------------------------------------------------------------------


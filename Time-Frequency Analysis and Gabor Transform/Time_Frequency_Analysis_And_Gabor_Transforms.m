close all;
clear all;

%Time analysis
%-------------
%Time domain
L = 10; %is the range of the time domain
%number of points in the time domain
n = 2048; 
%Define the time variable
t2 = linspace(0,L,n+1);
t = t2(1:n);
%Define signal
%arbitirary function/signal
S = (3*sin(2*t)+0.5*tanh(0.5*(t-3))+0.2*exp(-(t-4).^2)+1.5*sin(5*t)+4*cos(3*(t-6).^2))/10+(t/20).^3;
%Quadratic Chirp signal
% f0 = 50;
% f1 = 250;
% S = cos(2*pi*t.*(f0+(f1-f0)*t.^2/(3*L^2)));
% %electric guitar solo (audio file)
% [y,fs] = audioread('BoothBabyMantell.wav');
% %sound(y,fs);
% S = y(1:n)';
%plotting the signal
figure(1)
subplot(2,1,1)
plot(t,S)
title('Signal in time domain')
xlabel('Time')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
hold on

%Frequency analysis
%-----------------
%Define wave number
k = (2*pi./L).*[0:n/2-1 -n/2:-1];
ks = fftshift(k); %Fourier shift
%Fourier Transform of the signal
St = fftshift(fft(S));
%plotting
subplot(2,1,2)
plot(ks,abs(St)/max(abs(St))); %normalised
title('Signal in frequency domain')
xlabel('Frequency')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)

%Time frequency analysis (spectrogram)
%----------------------------------------
%A spectrogram is a visual representation of the spectrum of frequencies of a signal as it varies with time. 
%Initialise spectrogram matrix
Sgt_spec = []; 
%Define the time variable for iteration
tslide = 0:0.1:L; 
%Loop over all time steps
for j=1:length(tslide)

%Gabor Transform/window function
%option1
%g = exp(-(t-tslide(j)).^2); 
%option2: boarder window funciton
%g = exp(-0.2*(t-tslide(j)).^2);
%option3: narrower window funciton
g = exp(-5*(t-tslide(j)).^2);

%The broader the window, we have more accurate on frequency content 
%at the sacrifice of less accurate localisation on where the signal is in
%time. So, we can trust the information along the Y domain rather than the
%informtion o the X domain.
%the narrower the window, we have more accurate information on the time
%domain (time localisation) and less information on the frequency.

%Gabor Signal
Sg = g.*S; 

%Fourier Transofrm of the resulting signal
Sgt = fft(Sg);

%stor the spectrogram data at each iteration
Sgt_spec = [Sgt_spec;abs(fftshift(Sgt))/max(abs(Sgt))]; %matrix of spectrogram
%with rows to the time snapshots and coloums of Gabor signal

%plotting
figure(2)
%the signal and the sliding window function
subplot(3,1,1)
plot(t,S,'k',t,g,'r')
title('Signal with Gabor function in time domain')
xlabel('Time')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
%the convolution of the signal with the window function
subplot(3,1,2)
plot(t,Sg,'k')
title('Multiplication of the Signal with Gabor Transform in time domain')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
ylim([-1 1])
%the normalised resulting signal in frequency domain
subplot(3,1,3)
plot(ks,abs(fftshift(Sgt))/max(abs(Sgt)))
%axis([-50 50 0 1])
title('Fourier Transofrm of the new signal in frequency domain')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)

drawnow

pause(0.1)
end

%plotting the spectrogram
figure(3)
pcolor(tslide,ks,Sgt_spec.') 
shading interp
set(gca, 'Ylim', [-100 100])
colormap(jet)
xlabel('Time')
ylabel('Frequency')
title('Amplitude of a particular frequency at a particular time')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)


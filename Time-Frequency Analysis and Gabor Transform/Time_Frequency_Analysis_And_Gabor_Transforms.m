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
%optione1: arbitirary function/signal
S = (3*sin(2*t)+0.5*tanh(0.5*(t-3))+0.2*exp(-(t-4).^2)+1.5*sin(5*t)+4*cos(3*(t-6).^2))/10+(t/20).^3;
%----------
%option2: Quadratic Chirp signal
% f0 = 50;
% f1 = 250;
% S = cos(2*pi*t.*(f0+(f1-f0)*t.^2/(3*L^2)));
%----------
%option3: electric guitar solo (audio file)
% [y,fs] = audioread('BoothBabyMantell.wav');
% %sound(y,fs);
% S = y(1:n)';
%Data pre-processing
% %centering the data
% dsdt = diff(S(:))./diff(t(:));
% %rescaling
% mindsdt  = min(dsdt(:));
% maxdsdt  = max(dsdt(:));
% meandsdt = mean(dsdt(:));
% stddsdt  = std(dsdt(:));
% %rescale data from 0 to 1
% %Da1 = (Da - minDa)/(maxDa - minDa);
% %rescale data from -1 to 1
% dsdt2 = (dsdt - meandsdt)./stddsdt;
% maxdsdt2 = max(dsdt2(:));
% dsdt2 = dsdt2/maxdsdt2;
% S = dsdt2';
% %Add new column
% S = [0 S];
%----------

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
St = fft(S);
%plotting
subplot(2,1,2)
plot(ks,abs(fftshift(St))/max(abs(fftshift(St)))); %normalised
title('Signal in frequency domain')
xlabel('Frequency (Hz)')
ylabel('Intensity')
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
width = 5;
g = exp(-width*(t-tslide(j)).^2);

%The broader the window/filter, we have more accurate on frequency content 
%at the sacrifice of less accurate localisation on where the signal is in
%time. So, we can trust the information along the Y domain rather than the
%informtion o the X domain. For broad window we get the low frequency component
%that narrow filter usually throws them out becasue low frequency components are 
%not locaslised in time.

%the narrower the window, we have more accurate information on the time
%domain (time localisation) and less information on the frequency, i.e., poor frequency
%resolution. In this kind of filters, a lot of frequency components are thrown away.

%Gabor Signal
Sg = g.*S; 

%Fourier Transofrm of the resulting signal
Sgt = fft(Sg);

%stor the spectrogram data at each iteration
Sgt_spec = [Sgt_spec;abs(fftshift(Sgt))/max(abs(fftshift(Sgt)))]; %matrix of spectrogram
%with rows to the time snapshots and coloums of Gabor signal

%plotting
figure(2)
%the signal and the sliding window function
subplot(3,1,1)
plot(t,S,'k',t,g,'r')
title('Signal with Gabor function in time domain')
xlabel('Time (sec)')
xlim([0 max(t)])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
%the convolution of the signal with the window function
subplot(3,1,2)
plot(t,Sg,'k')
xlabel('Time (sec)')
xlim([0 max(t)])
title('Convolution of the Signal with Gabor Transform in time domain')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
%the normalised resulting signal in frequency domain
subplot(3,1,3)
plot(ks,abs(fftshift(Sgt))/max(abs(fftshift(Sgt))))
xlabel('Frequency (Hz)')
title('Fourier Transofrm of the new signal in frequency domain')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)

drawnow

pause(0.1)
end

%plot the spectrogram along with its data
figure;
subplot(2,1,1)
plot(t,S,'b')
xlim([0 L])
xlabel('Time (sec)')
ylabel('Signal')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
subplot(2,1,2)
pcolor(tslide,ks,Sgt_spec')
shading interp
colormap jet
xlim([0 L])
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)

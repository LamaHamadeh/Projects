
close all;
clear all;

%space
L = 10;
%n = 2048;
n = 60000;
%time
t2 = linspace(0,L,n+1);
t = t2(1:n);

%arbitirary function/signal
%S = (3*sin(2*t)+0.5*tanh(0.5*(t-3))+0.2*exp(-(t-4).^2)+1.5*sin(5*t)+4*cos(3*(t-6).^2))/10+(t/20).^3;

%Quadratic Chirp signal
%f0 = 50;
%f1 = 250;
%S = cos(2*pi*t.*(f0+(f1-f0)*t.^2/(3*L^2)));

%electric guitar solo (audio file)
[y,fs] = audioread('dspafsx_mono.wav');
sound(y,fs);
S = y(1:n)';

%distinguish between hot and cold water
%distinguish between male and female voices
%matlab auio samples to try: https://uk.mathworks.com/help/audio/ug/sample-audio-files.html

%plotting
figure(1)
subplot(2,1,1)
plot(t,S)
title('Signal in time domain')
xlabel('Time')
hold on

%Time frequency analysis (spectrogram)
%A spectrogram is a visual representation of the spectrum of frequencies of a signal as it varies with time. 
%--------
%Wave number
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);
%Fourier Transform of the signal
St = fft(S);
%plotting
figure(1)
subplot(2,1,2)
plot(ks,abs(fftshift(St))); %no time information with all spectral information
title('Signal in frequency domain')
xlabel('Frequency')

%Gabor Transform
%window function
g = exp((t-5).^2);
Sg = g.*S; %Gabor Signal
Sgt = fft(Sg); %Fourier Transofrm of the resulting signal
%plotting
figure(2)
subplot(3,1,1)
plot(t,S,'k',t,g,'r')
title('Signal with Gabor function in time domain')
subplot(3,1,2)
plot(t,Sg,'k')
title('Multiplication of the Signal with Gabor Transform')
subplot(3,1,3)
plot(ks,abs(fftshift(Sgt))/max(abs(Sgt)))
axis([-50 50 0 1])
title('Fourier Transofrm of the new signal')

%sliding the window function along the signal
Sgt_spec = [];
tslide = 0:0.1:10; %increments in time
for j=1:length(tslide)
%Gabor Transform
%window function
%g = exp(-(t-tslide(j)).^2); 
%broader window function
%g = exp(-0.2*(t-tslide(j)).^2);
%The broader the window, we have more accurate on frequency content 
%at the sacrifice of less accurate localisation on where the signal is in
%time. So, we can trust the information along the Y domain rather than the
%informtion o the X domain.
%narrower window funciton
g = exp(-5*(t-tslide(j)).^2);
%the narrower the window, we have more accurate information on the time
%domain (time localisation) and less information on the frequency.

Sg = g.*S; %Gabor Signal
Sgt = fft(Sg); %Fourier Transofrm of the resulting signal

Sgt_spec = [Sgt_spec;abs(fftshift(Sgt))]; %matrix of spectrogram
%with rows to the time snapshots and coloums of Gabor signal

%plotting
figure(3)
subplot(3,1,1)
plot(t,S,'k',t,g,'r')
title('Signal with Gabor function in time domain')
subplot(3,1,2)
plot(t,Sg,'k')
title('Multiplication of the Signal with Gabor Transform in time domain')
subplot(3,1,3)
plot(ks,abs(fftshift(Sgt))/max(abs(Sgt)))
axis([-50 50 0 1])
title('Fourier Transofrm of the new signal in frequency domain')
drawnow
pause(0.1)
end

%plotting the spectrogram
figure(4)
pcolor(tslide,ks,Sgt_spec.') 
shading interp
set(gca, 'Ylim', [-10000 10000])
%set(gca, 'Ylim', [-50 50])
colormap(jet)
xlabel('Time')
ylabel('Frequency')
title('Amplitude of a particular frequency at a particular time')




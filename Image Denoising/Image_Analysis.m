

close all;
clear all;

%read image
A = imread('Photo','jpeg'); %Glaxo Smith Kline Carbon Neutral Laboratory for Sustainable Chemistry
%plotting
subplot(2,2,1)
image(A)
axis square
set(gca, 'Xtick',[],'Ytick',[])
%convert the colour image from image to a matrix
A2 = double(A);
%define white noise to the whole image
noise = randn(3024, 3024,3);
%add the noise to our colour image
u = A2+50*noise;
%plotting
subplot(2,2,2)
image(uint8(u)) %show as an image
axis square
set(gca, 'Xtick',[],'Ytick',[])
%--------------------
%convert image to black and white
Abw = rgb2gray(A);
%plotting
subplot(2,2,3)
axis square
imshow(Abw);
colormap(gray(256));
%convert the BW image from image to a matrix
Abw= double(Abw);
%define white noise to the BW image
noise2 = randn(3024,3024);
%add the noise to our BW image
u2 = Abw+50*noise2;
%plotting
subplot(2,2,4)
image(uint8(u2)) %show as an image
axis square
set(gca, 'Xtick',[],'Ytick',[])
%--------------------
%Fourier Transform
%flipping the BW photo upside down
%Abw2 = Abw(3024:-1:1,:); 
Abw2 = flipud(Abw);
%plotting
figure(2)
subplot(1,2,1)
pcolor(Abw2)
axis square
shading interp
colormap(hot)
set(gca, 'Xtick',[],'Ytick',[])
%Apply Fourier Transform on the flipped BW photo
Abwt = abs(fftshift(fft2(Abw2)));
%plotting
subplot(1,2,2)
pcolor(log(Abwt))
axis square
shading interp
colormap(hot)
set(gca, 'Xtick',[],'Ytick',[])
%--------------------








close all;
clear all;

%read image
A = imread('Photo','jpeg'); %Glaxo Smith Kline Carbon Neutral Laboratory for Sustainable Chemistry
%convert the colour image from image format to a matrix
A2 = double(A);
%--------------------

%convert image to black and white
Abw = rgb2gray(A);
%convert the BW image from image format to a matrix
Abw= double(Abw);
%--------------------

%add noise to the BW image
B = Abw+100*randn(size(Abw,1),size(Abw,1));
%plotting
subplot(2,2,1)
imshow(uint8(B))
colormap(gray(256))
title('Noisy image')
%--------------------

%Filtering

%Fourier Transform
Bt = fft2(B);
%Shifting
Bts = fftshift(Bt);
%plotting
subplot(2,2,2)
pcolor(log(abs(Bts)))
shading interp
colormap(gray(256))

%Define the wave number in 2D
kx = 1:size(Abw,1);
ky = 1:size(Abw,1);
[KX,KY] = meshgrid(kx,ky);

%build a filter
width = 50; %the smaller the width the better the filtering
%Guassian filter
%F = exp(-width*(KX-1513).^2-width*(KY-1513).^2); 
%Shannon filter
F = zeros(size(Abw,1),size(Abw,1)); 
Fx = 1513-width:1:1513+width; %Or: (size(Abw,1)/2)+1-width:1:(size(Abw,1)/2)+1+width
Fy = 1513-width:1:1513+width;
F(Fx,Fy) = ones(length(2*width+1),length(2*width+1));
%apply/convolving the filter on/with the transform shifted image
Btsf = Bts.*F; 
%plotting
subplot(2,2,3)
pcolor(log(abs(Btsf)))
shading interp
colormap(gray(256));
%shift it back
Btf = ifftshift(Btsf);
%apply the inverse of Frourier Transform
Bf = ifft2(Btf);
%plotting
subplot(2,2,4)
imshow(uint8(Bf)) %show it as an image
colormap(gray(256));
title('Reconstructed Filtered Image')
%--------------------
%show the performance of the Gaussian filter as the width decreases
% figure(2)
% 
% width =[0.01 0.001 0.0001 0]; %For Gaussian filter
% 
% for j =1:length(width)
%     F = exp(-width(j)*(KX-1513).^2-width(j)*(KY-1513).^2);
%     Btsf = Bts.*F; %transformed, shifted, filtered
%     Btf = ifftshift(Btsf); %shift it back
%     Bf = ifft2(Btf); %bring it back to time domain
%     imshow(uint8(Bf))
%     colormap(gray(256))
%     drawnow
% end


%show the performance of the Shannon filter as the width decreases
% figure(2)
% 
% width = [10 50 100 200]; %For Shannon filter
% 
% for j =1:length(width)
%     F = zeros(3024,3024); 
%     F(1513-width(j):1:1513+width(j),1513-width(j):1:1513+width(j))...
%             = ones(length(2*width(j)+1),length(2*width(j)+1));
%     Btsf = Bts.*F; %transformed, shifted, filtered
%     Btf = ifftshift(Btsf); %shift it back
%     Bf = ifft2(Btf); %bring it back to time domain
%     imshow(uint8(Bf))
%     colormap(gray(256))
%     drawnow
% end


%--------------------


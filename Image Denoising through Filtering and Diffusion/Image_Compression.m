

close all;
clear all;

%read image
%A = imread('GSK','jpeg'); %Glaxo Smith Kline Carbon Neutral Laboratory for Sustainable Chemistry (UoN)
A = imread('Arkwright','png'); %Arkwright Building (NTU)
%plotting
subplot(2,3,1)
image(A)
axis square
set(gca, 'Xtick',[],'Ytick',[])
title('Original colour image')

%convert image to black and white
Abw = rgb2gray(A);
%plotting
subplot(2,3,2)
imshow(Abw);
axis square
title('Original BW image')

%Fourier Transform on the BW photo
Abw= double(Abw);
At = fft2(Abw);
Abwt = abs(fftshift(At));
%plotting
subplot(2,3,3)
pcolor(log(Abwt))
axis square
shading interp
colormap(hot)
set(gca, 'Xtick',[],'Ytick',[])
title('Foruier Transform')

%Apply image compression: look for  the small frequencies and zeroing them out
[nx,ny] = size(Abw);
count_pic = 4;
for threshold = 0.1*[0.001 0.005 0.01] * max(max(abs(At)))
    ind = abs(At) > threshold;
    count = nx*ny - sum(sum(ind));
    Atlow = At.*ind;
    percent = 100 - count/(nx*ny)*100;
    Alow = uint8(ifft2(Atlow));
    figure(1)
    subplot(2,3,count_pic)
    imshow(Alow)
    axis square
    count_pic = count_pic+1;
    drawnow
    title([num2str(percent) '% of FFT basis'])
end

%plot the 2D plot of th image frequencies
figure(2)
%original image
subplot(1,2,1)
%Anew = imresize(Abw,.1);
surf(fliplr(Abw));
shading interp
axis square
title('Uncompressed Image')
%compressed image
subplot(1,2,2)
%Anewlow = imresize(Alow,.1);
surf(fliplr(Alow));
shading interp
axis square
title('Compressed Image')


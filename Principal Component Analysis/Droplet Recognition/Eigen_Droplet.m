clear all;
close all;

%Training

%Read images
%Read images from Time 1
A1 = imread('Training\Time 1\1.tif'); %471   471     3
A1 = imresize(double(A1(:,:,1)),[60 60]); %choose only the red channel from 0 to 256
A2 = imread('Training\Time 1\2.tif'); %467   467     3
A2 = imresize(double(A2(:,:,1)),[60 60]);
A3 = imread('Training\Time 1\3.tif'); %477   477     3
A3 = imresize(double(A3(:,:,1)),[60 60]);
A4 = imread('Training\Time 1\4.tif'); %519   459     3
A4 = imresize(double(A4(:,:,1)),[60 60]);
%Read images from Time 2
B1 = imread('Training\Time 2\1.tif'); %471   471     3
B1 = imresize(double(B1(:,:,1)),[60 60]);
B2 = imread('Training\Time 2\2.tif'); %467   467     3
B2 = imresize(double(B2(:,:,1)),[60 60]);
B3 = imread('Training\Time 2\3.tif'); %477   477     3
B3 = imresize(double(B3(:,:,1)),[60 60]);
B4 = imread('Training\Time 2\4.tif'); %519   459     3
B4 = imresize(double(B4(:,:,1)),[60 60]);
%Read images from Time 3
C1 = imread('Training\Time 3\1.tif'); %471   471     3
C1 = imresize(double(C1(:,:,1)),[60 60]);
C2 = imread('Training\Time 3\2.tif'); %467   467     3
C2 = imresize(double(C2(:,:,1)),[60 60]);
C3 = imread('Training\Time 3\3.tif'); %477   477     3
C3 = imresize(double(C3(:,:,1)),[60 60]);
C4 = imread('Training\Time 3\4.tif'); %519   459     3
C4 = imresize(double(C4(:,:,1)),[60 60]);
%Read images from Time 4
D1 = imread('Training\Time 4\1.tif'); %471   471     3
D1 = imresize(double(D1(:,:,1)),[60 60]);
D2 = imread('Training\Time 4\2.tif'); %467   467     3
D2 = imresize(double(D2(:,:,1)),[60 60]);
D3 = imread('Training\Time 4\3.tif'); %477   477     3
D3 = imresize(double(D3(:,:,1)),[60 60]);
D4 = imread('Training\Time 4\4.tif'); %519   459     3
D4 = imresize(double(D4(:,:,1)),[60 60]);

%show images
figure (1)
subplot(4,4,1)
pcolor(flipud(A1)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,2)
pcolor(flipud(A2)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,3)
pcolor(flipud(A3)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,4)
pcolor(flipud(A4)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,5)
pcolor(flipud(B1)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,6)
pcolor(flipud(B2)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,7)
pcolor(flipud(B3)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,8)
pcolor(flipud(B4)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,9)
pcolor(flipud(C1)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,10)
pcolor(flipud(C2)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,11)
pcolor(flipud(C3)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,12)
pcolor(flipud(C4)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,13)
pcolor(flipud(D1)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,14)
pcolor(flipud(D2)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,15)
pcolor(flipud(D3)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,16)
pcolor(flipud(D4)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])

%Average droplets
AveA = (A1+A2+A3+A4)/4;
AveB = (B1+B2+B3+B4)/4;
AveC = (C1+C2+C3+C4)/4;
AveD = (D1+D2+D3+D4)/4;

%show the average
figure (2)
subplot(2,2,1)
pcolor(flipud(AveA)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(2,2,2)
pcolor(flipud(AveB)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(2,2,3)
pcolor(flipud(AveC)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(2,2,4)
pcolor(flipud(AveD)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])

%Data matrix
DB = [reshape(A1,1,60*60)
    reshape(A2,1,60*60)
    reshape(A3,1,60*60)
    reshape(A4,1,60*60)
    reshape(B1,1,60*60)
    reshape(B2,1,60*60)
    reshape(B3,1,60*60)
    reshape(B4,1,60*60)
    reshape(C1,1,60*60)
    reshape(C2,1,60*60)
    reshape(C3,1,60*60)
    reshape(C4,1,60*60)
    reshape(D1,1,60*60)
    reshape(D2,1,60*60)
    reshape(D3,1,60*60)
    reshape(D4,1,60*60)];

%construct covariance matrix
Cov = (DB')*(DB); 

%diagnolise covariance matrix, i.e., compute eigenvalues and eigenvectors
[V,E] = eigs(Cov,20,'lm');

%show the eigenvectors and the eigenvalues
figure(3)
subplot(2,3,1), drop1 = reshape(V(:,1), 60,60); pcolor(flipud(drop1)), shading interp, colormap gray 
% V(:,1) the first column of matrix V: the first eigenvector which contains the main characteristics all images have.
subplot(2,3,2), drop2 = reshape(V(:,2), 60,60); pcolor(flipud(drop2)), shading interp, colormap gray
subplot(2,3,3), drop3 = reshape(V(:,3), 60,60); pcolor(flipud(drop3)), shading interp, colormap gray
subplot(2,3,4), drop4 = reshape(V(:,4), 60,60); pcolor(flipud(drop4)), shading interp, colormap gray
subplot(2,3,5), drop5 = reshape(V(:,5), 60,60); pcolor(flipud(drop5)), shading interp, colormap gray
subplot(2,3,6), plot(diag(E),'ko', 'Linewidth',[2]) %try semilogx on the y axis
set(gca, 'Fontsize', [14])

%transform the average images from matrices to vectors
vecA = reshape(AveA,1,60*60);
vecB = reshape(AveB,1,60*60);
vecC = reshape(AveC,1,60*60);
vecD = reshape(AveD,1,60*60);

%project the average images on the eigenvectors, i.e., classification functions.
projA = vecA * V;
projB = vecB * V;
projC = vecC * V;
projD = vecD * V;

%show the projection of the average images on the 20 eigenvectors
figure(4)
subplot(2,2,1), bar(projA(1:20)), set(gca, 'Xlim',[0 20], 'Ylim', [-2000 2000], 'Xtick', [], 'Ytick', [])
text(12,-1700, 'Time 1', 'Fontsize', [15])
subplot(2,2,2), bar(projB(1:20)), set(gca, 'Xlim',[0 20], 'Ylim', [-2000 2000], 'Xtick', [], 'Ytick', [])
text(12,-1700, 'Time 2', 'Fontsize', [15])
subplot(2,2,3), bar(projC(1:20)), set(gca, 'Xlim',[0 20], 'Ylim', [-2000 2000], 'Xtick', [], 'Ytick', [])
text(12,-1700, 'Time 3', 'Fontsize', [15])
subplot(2,2,4), bar(projD(1:20)), set(gca, 'Xlim',[0 20], 'Ylim', [-2000 2000], 'Xtick', [], 'Ytick', [])
text(12,-1700, 'Time 4', 'Fontsize', [15])

%-----------------------------------

%Testing

%Read images
T1 = imread('Test\Time 1\5.tif');
T1 = imresize(double(T1(:,:,1)),[60 60]);
T2 = imread('Test\Time 2\5.tif');
T2 = imresize(double(T2(:,:,1)),[60 60]);
T3 = imread('Test\Time 3\5.tif');
T3 = imresize(double(T3(:,:,1)),[60 60]);
T4 = imread('Test\Time 4\5.tif');
T4 = imresize(double(T4(:,:,1)),[60 60]);

%reshape them from matrices to vectors
vecT1 = reshape(T1,1,60*60);
vecT2 = reshape(T2,1,60*60);
vecT3 = reshape(T3,1,60*60);
vecT4 = reshape(T4,1,60*60);

%project the test data(vectors) on the eigenvectors space
projT1 = vecT1*V;
projT2 = vecT2*V;
projT3 = vecT3*V;
projT4 = vecT4*V;

%reconstruct the projeccted data back into an image
recon1 = V*projT1'; rec1 = reshape(recon1, 60 ,60);
recon2 = V*projT2'; rec2 = reshape(recon2, 60 ,60);
recon3 = V*projT3'; rec3 = reshape(recon3, 60 ,60);
recon4 = V*projT4'; rec4 = reshape(recon4, 60 ,60);

%show prediction results
figure(5)
subplot(4,4,1)
pcolor(flipud(T1)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,2)
bar(projT1(1:20)), set(gca, 'Xlim',[0 20], 'Ylim', [-2000 2000], 'Xtick', [], 'Ytick', [])
subplot(4,4,3)
pcolor(flipud(rec1)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,4)
E1 = [norm(projA-projT1)/norm(projA) norm(projB-projT1)/norm(projB) norm(projC-projT1)/norm(projC) norm(projD-projT1)/norm(projD)];
bar(E1,'b','EdgeColor','none') 

subplot(4,4,5)
pcolor(flipud(T2)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,6)
bar(projT2(1:20)), set(gca, 'Xlim',[0 20], 'Ylim', [-2000 2000], 'Xtick', [], 'Ytick', [])
subplot(4,4,7)
pcolor(flipud(rec2)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,8)
E2 = [norm(projA-projT2)/norm(projA) norm(projB-projT2)/norm(projB) norm(projC-projT2)/norm(projC) norm(projD-projT2)/norm(projD)];
bar(E2,'b','EdgeColor','none') 

subplot(4,4,9)
pcolor(flipud(T3)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,10)
bar(projT3(1:20)), set(gca, 'Xlim',[0 20], 'Ylim', [-2000 2000], 'Xtick', [], 'Ytick', [])
subplot(4,4,11)
pcolor(flipud(rec3)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,12)
E3 = [norm(projA-projT3)/norm(projA) norm(projB-projT3)/norm(projB) norm(projC-projT3)/norm(projC) norm(projD-projT3)/norm(projD)];
bar(E3,'b','EdgeColor','none') 

subplot(4,4,13)
pcolor(flipud(T4)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,14)
bar(projT4(1:20)), set(gca, 'Xlim',[0 20], 'Ylim', [-2000 2000], 'Xtick', [], 'Ytick', [])
subplot(4,4,15)
pcolor(flipud(rec4)), shading interp, colormap(gray), set(gca, 'Xtick',[],'Ytick',[])
subplot(4,4,16)
E4 = [norm(projA-projT4)/norm(projA) norm(projB-projT4)/norm(projB) norm(projC-projT4)/norm(projC) norm(projD-projT4)/norm(projD)];
bar(E4,'b','EdgeColor','none') 
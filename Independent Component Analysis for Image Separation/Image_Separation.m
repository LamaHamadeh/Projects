
close all;
clear all;


%Read images taken from my iphone 
S1 = imread('IS1','jpeg'); 
S2 = imread('IS2','jpeg'); 

%Show raw images
figure
subplot(3,2,1), imshow(S1), title('Centre for Sustainable Chemistry')
subplot(3,2,2), imshow(S2), title('My frontdoor wreath')
hold on

%Mix images twice to have two sets of measurements
A = [4/5 3/5 ; 1/2 2/3]; %arbitarty mixing matrix
% A = [4/5 3/5 ; 1/2 2/3];

x1 = A(1,1)*S1+A(1,2)*S2; %first measurement 
x2 = A(2,1)*S1+A(2,2)*S2; %second measurement

%transfer the mixing images to double precision
x1 = double(x1);
x2 = double(x2);

%show miximg measurements/images
subplot(3,2,3), imshow(uint8(x1))
subplot(3,2,4), imshow(uint8(x2))


%SVD Reconstruction process
%--------------------------
%check the size of S1 and S2
[m1,n1] = size(x1); %3024 9072
[m2,n2] = size(x2); %3024 9072
%reshape of mixing images to column vectors
x1 = reshape(x1, m1*n1, 1);
x2 = reshape(x2, m2*n2, 1);

%step 1: maximal/minimal variance angle detection
%We assume a mean-zero distribution
x1 = x1-mean(x1);
x2 = x2-mean(x2);
%critical angle that refers to the maximum/minimum variance
theta0 = 0.5*atan(-2*sum(x1.*x2)/sum(x1.^2-x2.^2));
%the unitary/rotation matrix
Us=[cos(theta0) sin(theta0) ; -sin(theta0) cos(theta0)];
% %-------------

%step 2: scaling of the principal components
%calculate variance in the principle direction
sig1 = sum((x1*cos(theta0)+x2*sin(theta0)).^2);
%%calculate variance in the secondary direction
sig2 = sum((x1*cos(theta0-pi/2)+x2*sin(theta0-pi/2)).^2);
%construct sigma matrix
sigma = [1/sqrt(sig1) 0 ; 0 1/sqrt(sig2)];
% %-------------

%step 3: Rotation to separability (make probability density separable)
%construct the images that have undergone the two steps of transformation
x1bar = sigma(1,1)*(Us(1,1)*x1+Us(1,2)*x2); %xbar = sigma^-1*Ustar*x => S=V.xbar
x2bar = sigma(2,2)*(Us(2,1)*x1+Us(2,2)*x2);

% x1bar = reshape(x1bar, m1*n1,1);
% x2bar = reshape(x2bar, m2*n2,1);

%critical angle for the third step
phi0 = 0.25*atan(-sum(2*(x1bar.^3).*x2bar-2*x1bar.*(x2bar.^3))...
    /sum(3*(x1bar.^2).*(x2bar.^2)-0.5*(x1bar.^4)-0.5*(x2bar.^4)));
%the rotation matrix 
V=[cos(phi0) sin(phi0); -sin(phi0) cos(phi0)];
%get the approximated final separated images
s1bar=V(1,1)*x1bar+V(1,2)*x2bar;
s2bar=V(2,1)*x1bar+V(2,2)*x2bar;
% %-------------
%plotting/showing the spearated images
%rescales the spearated images to values ranging from 0 to 255 
%first image
min1=min(min(min(s1bar)));
s1bar=s1bar-min1;
max1=max(max(max(s1bar)));
s1bar=s1bar*(255/max1);
%second image
min2=min(min(min(s2bar)));
s2bar=s2bar-min2;
max2=max(max(max(s2bar)));
s2bar=s2bar*(255/max2);
%reshape the vectors back to images
s1bar = reshape(s1bar, 3024, 3024, 3);
s2bar = reshape(s2bar, 3024, 3024, 3);
%show the approximated separated images
subplot(3,2,5), imshow(uint8(s1bar))
subplot(3,2,6), imshow(uint8(s2bar))
% %-------------

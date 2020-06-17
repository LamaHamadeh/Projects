

clear all;
close all;

%read image
A = imread('Durham_Castle.jpg');

%convertt to gray
X = double(rgb2gray(A));
%define image dimensions
nx = size(A,1);
ny = size(A,2);

%show original image
figure;
subplot(2,2,1)
imagesc(X)
axis off
axis square
colormap gray
title('Original')

%Apply SVD
[U,S,V] = svd(X);

%truncate SVD components: compress image
plotind = 2;
for r = [5 20 100] %truncation value
    Xapprox = U(:,1:r)*S(1:r,1:r)*V(:,1:r)'; %approximated image
    %show the approximated image
    subplot(2,2,plotind)
    plotind = plotind + 1;
    imagesc(Xapprox)
    axis off
    axis square
    title(['r= ',num2str(r)]);
end

%to show eigenvalues entire spectrum
figure
subplot(1,2,1)
semilogy(diag(S),'k','LineWidth',2)
xlabel('Number of Eigenvalues, r')
ylabel('Magnitude of Eigenvalues')
grid on
axis square
xlim([-50 1550])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
%Cumulative Energy
subplot(1,2,2)
plot(cumsum(diag(S))/sum(diag(S)),'k','LineWidth',2)
xlabel('Number of Eigenvalues, r')
ylabel('Cumulative Energy')
grid on
axis square
xlim([-50 1550])
ylim([0 1.1])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
% the first column in U and V accounts for 28% of the total imformation in 
% the dog.
%------------------------------------





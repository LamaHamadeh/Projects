
%Case 1: this is the ideal case where you consider a small displacement of
% the mass in the z direction and the ensuing oscillations. In this case, 
% the entire motion is in the z direction with the simple harmonic motion 
% being observed (camN_1.mat where N=1,2,3).

close all;
clear all;

%load the datasets
%first dataset from the first camera
ds1 = load('cam1_2.mat');
%second dataset from the second camera
ds2 = load('cam2_2.mat');
%third dataset from the third camera
ds3 = load('cam3_2.mat');
%------------------

%check the size of each dataset
[row1,col1,rgb1,frame1] = size(ds1.vidFrames1_2);
[row2,col2,rgb2,frame2] = size(ds2.vidFrames2_2);
[row3,col3,rgb3,frame3] = size(ds3.vidFrames3_2);
%each dataset is a 4D matrix. The first two dimensions represent the number
%of rows and columns of the images captured from the recorded videos.
%The third dimension represents the three colour channels: RGB. The fourth
%dimension represents the number of images/frames in that dataset.
%------------------

%Image Processing & Can Movement Tracking 
%Camera 1
%---------
%define the limits where we can see (by eye) the can lies within it 
left_edge1   = 200;
right_edge1  = 450;
top_edge1    = 450; 
bottom_edge1 = 150;

%initialise x and y location matrices
X1_2 = zeros(1,frame1); % Store object's x location
Y1_2 = zeros(1,frame1); % Store object's y location

for i = 1:frame1
    % Chnage each frame's type to a "double" instead of its uint8
    % format
    img1(:,:,:,i) = double(rgb2gray(ds1.vidFrames1_2(:,:,:,i)));
    % Image Cropping: isolate the paint can in the video by making 
    % areas around the can black, i.e., pixels values equal to zero.
    for r = 1:row1
        for c = 1:col1
            if (r > top_edge1) || (r < bottom_edge1) ...
                    || (c > right_edge1) || (c < left_edge1)
                img1(r,c,i) = 0;
            end

        end
    end
    % Object Tracking
    %find pixel with maximum itensity that corresponds to light
    %flash of the camera and its pixel index
    [maximg1, indx1] = max(max(img1(:,:,i)));
    %cFind out the row and column of the indices of all maximum points.
    %Thsi corresponds to the x and y cooridnates of the maximum point
    [maxx1,maxy1] = ind2sub([row1 col1],indx1);
    %record all x-y locations of maximum points
    X1_2(i) = maxx1;
    Y1_2(i) = maxy1;
end

% %show the image
% k = 15;
% for j = 1:k
%     pcolor(img1(:,:,j))
%     shading interp
%     colormap("gray")
%     yline(top_edge1,'r','LineWidth',2)
%     yline(bottom_edge1,'r','LineWidth',2)
%     xline(left_edge1,'b','LineWidth',2)
%     xline(right_edge1,'b','LineWidth',2)
%     xlabel('Width','FontSize',14)
%     ylabel('Height','FontSize',14)
%     set(gca,'TickLabelInterpreter','latex')
%     set(gca,'YDir','normal')
%     hold on
%     pause(0.1)
% end
%------------------

%Camera 2
%---------
%define the limits where we can see (by eye) the can lies within it 
left_edge2   = 200;
right_edge2  = 450;
top_edge2    = 450; 
bottom_edge2 = 150;

%initialise x and y location matrices
X2_2 = zeros(1,frame2); % Store object's x location
Y2_2 = zeros(1,frame2); % Store object's y location

for i = 1:frame2
    % Chnage each frame's type to a "double" instead of its uint8
    % format
    img2(:,:,:,i) = double(rgb2gray(ds2.vidFrames2_2(:,:,:,i)));
    % Image Cropping: isolate the paint can in the video by making 
    % areas around the can black, i.e., pixels values equal to zero.
    for r = 1:row2
        for c = 1:col2
            if (r > top_edge2) || (r < bottom_edge2) ...
                    || (c > right_edge2) || (c < left_edge2)
                img2(r,c,i) = 0;
            end

        end
    end
    % Object Tracking
    %find pixel with maximum itensity that corresponds to light
    %flash of the camera and its pixel index
    [maximg2, indx2] = max(max(img2(:,:,i)));
    %cFind out the row and column of the indices of all maximum points.
    %Thsi corresponds to the x and y cooridnates of the maximum point
    [maxx2,maxy2] = ind2sub([row2 col2],indx2);
    %record all x-y locations of maximum points
    X2_2(i) = maxx2;
    Y2_2(i) = maxy2;
end

% %show the image
% k = 5;
% for j = 1:k
%     pcolor(img2(:,:,j))
%     shading interp
%     colormap("gray")
%     yline(top_edge2,'r','LineWidth',2)
%     yline(bottom_edge2,'r','LineWidth',2)
%     xline(left_edge2,'b','LineWidth',2)
%     xline(right_edge2,'b','LineWidth',2)
%     xlabel('Width','FontSize',14)
%     ylabel('Height','FontSize',14)
%     set(gca,'TickLabelInterpreter','latex')
%     set(gca,'YDir','normal')
%     hold on
%     pause(0.1)
% end
%------------------

%Camera 3
%---------
%define the limits where we can see (by eye) the can lies within it 
left_edge3   = 200;
right_edge3  = 450;
top_edge3    = 450; 
bottom_edge3 = 150;

%initialise x and y location matrices
X3_2 = zeros(1,frame3); % Store object's x location
Y3_2 = zeros(1,frame3); % Store object's y location

for i = 1:frame3
    % Chnage each frame's type to a "double" instead of its uint8
    % format
    img3(:,:,:,i) = double(rgb2gray(ds3.vidFrames3_2(:,:,:,i)));
    % Image Cropping: isolate the paint can in the video by making 
    % areas around the can black, i.e., pixels values equal to zero.
    for r = 1:row3
        for c = 1:col3
            if (r > top_edge3) || (r < bottom_edge3) ...
                    || (c > right_edge3) || (c < left_edge3)
                img3(r,c,i) = 0;
            end

        end
    end
    % Object Tracking
    %find pixel with maximum itensity that corresponds to light
    %flash of the camera and extract its pixel index
    [maximg3, indx3] = max(max(img3(:,:,i)));
    %Find out the row and column of the indices of all maximum points.
    %This corresponds to the x and y cooridnates of the maximum points
    [maxx3,maxy3] = ind2sub([row3 col3],indx3);
    %record all x-y locations of maximum points
    X3_2(i) = maxx3;
    Y3_2(i) = maxy3;
end

% %show the image
% k = 15;
% for j = 1:k
%     pcolor(img3(:,:,j))
%     shading interp
%     colormap("gray")
%     yline(top_edge3,'r','LineWidth',2)
%     yline(bottom_edge3,'r','LineWidth',2)
%     xline(left_edge3,'b','LineWidth',2)
%     xline(right_edge3,'b','LineWidth',2)
%     xlabel('Width','FontSize',14)
%     ylabel('Height','FontSize',14)
%     set(gca,'TickLabelInterpreter','latex')
%     set(gca,'YDir','normal')
%     hold on
%     pause(0.1)
% end

%------------------

%Construct the experiement dataset
%data matrix 6x226
X = [X1_2          ; Y1_2           ; 
    X2_2(1:frame1) ; Y2_2(1:frame1) ;
    X3_2(1:frame1) ; Y3_2(1:frame1) ];
%extract the size of the dataset
[Xm,Xn] = size(X);
%------------------

%Principal Component analysis (PCA)
% Compute the mean of all samples
mn = mean(X,2);

% Build the the mean-centring dataset (Data-mean) by subtracting the
% %mean mn from the samples and rotate it.
% X = X-repmat(mn,1,Xn);
X = X - mn*ones(1,Xn);

% Construct the covariance matrix L:
L = X*X'; 
% L = cov(X);

% Calculate the eigenvectors and the eigenvalues of the covariance
% matrix L
[V,D] = eig(L); %Here: V are the eigenvectors and the diagonal elements of 
% D are the eigenvalues of the matrix L.
%extract the principal components eigenvalues, i.e.m diagonal elements of D
lambda1 = diag(D);

%sort in decreasing order
[dummy,m_arrange] = sort(-1*lambda1);
lambda1 = lambda1(m_arrange);
V = V(:,m_arrange);

% Get PCA ranks/scores/modes by projecting our mean_centred dataset
% on the covariance matrix eigenvectors (principal component vector space)
Y1 = V'*X;
%visulaise the scores/Modes
figure;
subplot(1,2,1)
plot(Y1(1,:),'b','LineWidth',1.5)
hold on
plot(Y1(2,:),'r','LineWidth',1.5)
hold on
plot(Y1(3,:),'g','LineWidth',1.5)
xlabel('Time/Frame','FontSize',14)
ylabel('Displacement Position','FontSize',14)
title('Displacement Time Evolution')
legend('Mode 1','Mode 2','Mode 3')
set(gca,'TickLabelInterpreter','latex')
axis square
grid on
hold on

%visualise the principal components
subplot(1,2,2)
plot(lambda1/sum(lambda1),'bo','linewidth',1.5)
xlabel('SVD/PCA Principal Components','FontSize',14)
ylabel('Variance in Percentage ','FontSize',14)
title('Principal Components of the Data')
set(gca,'TickLabelInterpreter','latex')
axis square
grid on
% ------------------

%Singular Value Decomposition
%calculate the eignvalues and both eigenvectors
[U,S,V] = svd(X/sqrt(Xn-1),'econ');
%extract the diagonals, i.e., eigenvalues
lambda2= diag(S).^2;
%calculate the SVD modes, i.e., projections
Y2 = U'*X;
%visualise the first three modes
figure;
subplot(1,2,1)
plot(Y2(1,:),'b','LineWidth',1.5)
hold on
plot(Y2(2,:),'r','LineWidth',1.5)
hold on
plot(Y2(3,:),'g','LineWidth',1.5)
xlabel('Time/Frame','FontSize',14)
ylabel('Displacement Position','FontSize',14)
title('Displacement Time Evolution')
legend('Mode 1','Mode 2','Mode 3')
set(gca,'TickLabelInterpreter','latex')
grid on
axis square
hold on
%visualise the principal components
subplot(1,2,2)
plot(lambda2/sum(lambda2),'bo','linewidth',1.5)
xlabel('SVD/PCA Principal Components','FontSize',14)
ylabel('Variance in Percentage ','FontSize',14)
title('Principal Components of the Data')
set(gca,'TickLabelInterpreter','latex')
grid on
axis square
%------------------

%Graph/dynamics Explanation:
% We see here that adding noise through camera shaking does lower the
% variance in the first principal component but overall as seen in the right
% graph, we can still see the primary oscillation from the first principal 
% component, even if it is noticeably more noisy. Clearly introdcuing noise
% ot the data has impacted the quality of the results and the performance
% of the PCA on reducing the dimensionality of the overal dyanmics.
%------------------


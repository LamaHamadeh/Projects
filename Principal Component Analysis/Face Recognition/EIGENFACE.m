clear all;
close all; 


% Size of each picture
m = 200; %number of rows
n = 175; %number of columns

% Number of sample pictures (per person)
N = 20;

avg = zeros(m*n,1);  % the average face with size mXn
A = []; %store everything in this matrix

%% Load Arnold Schwarzenegger

%show the images of Arnold one by one
count = 0;
for j = 1:N
    figure(1) 
    %file path
    ff = ['faces/Arnold',num2str(j,'%02d'),'.jpg'];
    %read image
    u = imread(ff); % Read the image into a matrix u
    %show image
    imshow(u)
    %one channle colour
    if(size(u,3)==1)
        M=double(u);
    else
        M=double(rgb2gray(u)); 
    end
    
    pause(0.1);
    %convert M (n,m) into a column vector (nXm,1)
    R = reshape(M,m*n,1);
    %store this vector column into the main matrix A
    A = [A,R];
    %compute average
   avg = avg + R;
   %go to the next image
   count = count + 1;
end

%% Load Sylvester Stallone

%show the images of Arnold one by one
count = 0;
for j = 1:N
    figure(1) 
    %file path
    ff = ['faces/Stallone',num2str(j,'%02d'),'.jpg'];
    %read image
    u = imread(ff); % Read the image into a matrix u
    %show image
    imshow(u)
    %one channle colour
    if(size(u,3)==1)
        M=double(u);
    else
        M=double(rgb2gray(u)); 
    end
    
    pause(0.1);
    %convert M (n,m) into a column vector (nXm,1)
    R = reshape(M,m*n,1);
    %store this vector column into the main matrix A
    A = [A,R];
    %compute average
   avg = avg + R;
   %go to the next image
   count = count + 1;
end

%% Calculate the "averaged" face

avgTS = uint8(reshape(avg,m,n)); %reshape the average column vector back 
%into a picture to show it
%show it
figure(1)
imshow(avgTS);


%% Center the sample pictures at the "origin"
%mean subtracting
figure(1)
for j = 1:2*N
    %subtract the average
    A(:,j) = A(:,j) - avg;
    %plot it
    R = reshape(A(:,j),m,n);
    imshow(R);

    pause(.1);
end

%in these images we notices the dominant features of both celebrities.
%without the details.

%%  Computing the SVD: the face print

[U,S,V] = svd(A,'econ'); %or [U,S,V] = svd(A,0);
%economy SVD is a lot faster in computing
Phi = U(:,1:2*N);
Phi(:,1) = -1*Phi(:,1);
figure(2)
count = 1;
for i=1:3
    for j=1:3
        subplot(3,3,count)
        imshow(uint8(25000*reshape(Phi(:,count),m,n)));
        count = count + 1;
    end
end


%% project each image onto basis 

for j = 1:N
    imvec = A(:,j);
    ARN(:,j) = imvec'*Phi(:,1:3); %project Arnold pictures on the basis
end
for j = 1:N
    imvec = A(:,N+j);
    STAL(:,j) = imvec'*Phi(:,1:3); %project Stalon pictures on the basis
end

figure(3)

plot3(ARN(1,:),ARN(2,:),ARN(3,:),'r.','MarkerSize',30)
hold on
plot3(STAL(1,:),STAL(2,:),STAL(3,:),'b.','MarkerSize',30)
xlabel('PCA 1') %eigenface 1
ylabel('PCA 2') %eigenface 2
zlabel('PCA 3') %eigenface 3
legend('STALLONE','TAYLOR')
%the clustering between both celebrities is clear and noticeable

%% add some unexpected pics: testing

u = imread('faces/teststallone1.jpg');        
figure(4)
subplot(1,2,1)
imshow(u);
u = double(rgb2gray(u));
ustal = reshape(u,m*n,1)-avg;
stalpts = ustal'*Phi(:,1:3); %project on the first 3 PCs
v = imread('faces/testterminator8.jpg');
subplot(1,2,2)
imshow(v);
v = double(rgb2gray(v));
vterm = reshape(v,m*n,1)-avg;
termpts = vterm'*Phi(:,1:3); %project on the first 3 PCs
%plotting
figure(3)
plot3(stalpts(1),stalpts(2),stalpts(3),'g.','MarkerSize',30)
plot3(termpts(1),termpts(2),termpts(3),'k.','MarkerSize',30)

% the terminator (black dot) is classified as a Arnold type of picture
% and the stalon (green dot) is clasified as a Stalon type picture.




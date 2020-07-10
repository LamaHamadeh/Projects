
%Eigenvalues and Eigenvectors/Eigenstates of the double-well with a Guassian barrier 
%using Finite Difference Method

clear all;
close all;

%Create the Hamiltonian matrix
dy = 0.01; %size of the step
ymin = -10; 
ymax = 10;
N = ((ymax-ymin)/dy)+1; %number of points
y = ymin:dy:ymax; %space variable
H = zeros(N,N); %initialise the Hamiltonian matrix

%double well parameters
v0 = 5; %the height (amplitude) of the Gaussian peak at the centre of the double well
b = 0.3; %the inverse of the width of the Gaussian peak at the centre of the double well

%diagonal elements
for n = 1:N %number of rows
    for m = 1:N %number of columns
        if (n==m)
            H(n,m) = (ymin+(n-1)*dy)^2/2+v0*exp(-b*(ymin+(n-1)*dy)^2)+1/(dy)^2; %the value of each diagonal element
        end
    end
end

%the off-diagonal elements
for n =1:N-1
    H(n,n+1) = -1/(2*(dy)^2); %the value of each upper off-diagonal element
end

for n = 2:N
    H(n,n-1) = -1/(2*(dy)^2); %the value of each lower off-diagonal element
end

%calculate the eignevalues and eigenstates
 [V,D] = eig(H);
 
%plotting the eigenstates
figure;
subplot(2,2,1)
plot(y,V(:,1),'b')
subplot(2,2,2)
plot(y,V(:,2),'r')
subplot(2,2,3)
plot(y,V(:,3),'k')
subplot(2,2,4)
plot(y,V(:,4),'g')

%checking the normalisation of the eignestate
Vprod = conj(V(:,3)).*(V(:,3))*dy; %a product between the eigenstate and 
%its conjugate
Vnorm = sum(Vprod); %the sum over all the region points
Vnormalised = V(:,3)/sqrt(Vnorm); %normalising the resulting eigenstate
 
%plotting
figure;
plot(y,Vnormalised) %the normalised eigenstate
%comparing between the numerical and analytical solution
Analytical_solution = (1/pi)^0.25*exp(-y.^2/2);
figure;
plot(y,Vnormalised, 'b',y,Analytical_solution,'ro')
 
%checking the othogonality of the ground state and the first eigenstate
%In this casse, the vanishing the product between these two eigenstates is
%not quite so obvious, but it follows from the fact that V(:,1) is even
%and V(:,2) an odd function of x, i.e., f(x) = -f(-x).
%Thus, their product is an odd fucntion of x, the integral of an odd
%function from-A to +A is zero.
Vprod = conj(V(:,2)).*(V(:,1))*dy;
Vorth = sum(Vprod);
 
%the probability of finding a particle in a certain energy level
pro = conj(V(:,50)).*(V(:,50))*dy;
figure;
plot(y,pro)
 
%looking at the eignevalues
e0 = eig(H);
figure;
plot(e0)



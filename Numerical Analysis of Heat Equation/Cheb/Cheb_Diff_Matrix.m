
%Matlab function to compute differentiation matrix D_N taken from the book:
%Spectral Methods in MATLAB
%By Lloyd N. Trefethen

%another source:
%Nathan Kutz YouTube channel
%https://www.youtube.com/watch?v=i95pnk7na1M&t=35s

%define a function that returns matrix D and x 
%and needs N: number of Fourier points
function [D,x] = Cheb_Diff_Matrix(N)

if N==0 , D=0; x=1; return, end
%define x
x = cos(pi*(0:N)/N)';
%define matrix coefficients
c = [2;ones(N-1,1);2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
%building D
D = (c*(1./c)')./(dX+(eye(N+1))); %off-diagonals
D = D-diag(sum(D,2)); %diagonals





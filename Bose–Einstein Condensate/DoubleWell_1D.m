
%Eigenvalues and Eigenvectors/Eigenstates of the double well potential trap 
%using Finite Difference Method
%Time independent Schrodinger Equation

clear all;
close all;

%Create the space vaiable
dy = 0.01; %size of the step
ymin = -10; %the minimum of the domain
ymax = 10; %the maximum of the domain
N = ((ymax-ymin)/dy)+1; %number of points in the domain
y = ymin:dy:ymax; %space variable
%-------------

%double-well potential parameters
%option1
% v0 = 9; %the height (amplitude) of the Gaussian peak at the centre of the double-well 
% % b = 0.25; %the inverse of the width of the Gaussian peak at the centre of the double-well.
% b = 1.5; 
% Vdw = y.^2/2+v0*exp(-b.*y.^2);
%option2
b = 0.5;
%gamma = 0.044877;
gamma = 0.044;
Vdw = -b*y.^2+gamma*y.^4+1.393;
%plotting
figure;
plot(y,Vdw,'k','LineWidth',2)
xlabel('$y$','Interpreter','latex')
ylabel('$V_{\mathrm{dw}(y)}$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
drawnow
%-------------

%Create the Hamiltonian matrix
H = zeros(N,N); %initialise the Hamiltonian matrix
%diagonal elements
for n = 1:N %number of rows
    for m = 1:N %number of columns
        if (n==m)
            H(n,m) = 2/dy^2+Vdw(n); %yn = ymin+(n-1)*dy
        end
    end
end
%the off-diagonal elements
for n =1:N-1
    H(n,n+1) = -1/dy^2; %the value of each upper off-diagonal element
end

for n = 2:N
    H(n,n-1) = -1/dy^2; %the value of each lower off-diagonal element
end
%-------------

%calculate the eignevalues and eigenstates
 [Eigen_Fun,Eigen_Val] = eig(H);
%-------------
                            %%%Eigenstates%%%

%plotting the eigenstates
figure;
%ground state
subplot(2,2,1)
%checking the normalisation of the eignestate
Vprod1 = conj(Eigen_Fun(:,1)).*(Eigen_Fun(:,1))*dy; %a product between the eigenstate and 
%its conjugate
Vnorm1 = sum(Vprod1); %the sum over all the region points
Vnormalised1 = Eigen_Fun(:,1)/sqrt(Vnorm1); %normalising the resulting eigenstate
plot(y,Vnormalised1,'k','LineWidth',2)
xlabel('$y$','Interpreter','latex')
ylabel('$\psi_0(y)$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
%---
%first excited state
subplot(2,2,2)
%checking the normalisation of the eignestate
Vprod2 = conj(Eigen_Fun(:,2)).*(Eigen_Fun(:,2))*dy; %a product between the eigenstate and 
%its conjugate
Vnorm2 = sum(Vprod2); %the sum over all the region points
Vnormalised2 = Eigen_Fun(:,2)/sqrt(Vnorm2); %normalising the resulting eigenstate
plot(y,Vnormalised2,'k','LineWidth',2)
xlabel('$y$','Interpreter','latex')
ylabel('$\psi_1(y)$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
%---
%second excited state
subplot(2,2,3)
%checking the normalisation of the eignestate
Vprod3 = conj(Eigen_Fun(:,3)).*(Eigen_Fun(:,3))*dy; %a product between the eigenstate and 
%its conjugate
Vnorm3 = sum(Vprod3); %the sum over all the region points
Vnormalised3 = Eigen_Fun(:,3)/sqrt(Vnorm3); %normalising the resulting eigenstate
plot(y,Vnormalised3,'k','LineWidth',2)
xlabel('$y$','Interpreter','latex')
ylabel('$\psi_2(y)$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
subplot(2,2,4)
%---
%third excited state
%checking the normalisation of the eignestate
Vprod4 = conj(Eigen_Fun(:,4)).*(Eigen_Fun(:,4))*dy; %a product between the eigenstate and 
%its conjugate
Vnorm4 = sum(Vprod4); %the sum over all the region points
Vnormalised4 = Eigen_Fun(:,4)/sqrt(Vnorm4); %normalising the resulting eigenstate
plot(y,Vnormalised4,'k','LineWidth',2)
xlabel('$y$','Interpreter','latex')
ylabel('$\psi_3(y)$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
%---

%checking the othogonality of the ground state and the first eigenstate
%In this casse, the vanishing the product between these two eigenstates is
%not quite so obvious, but it follows from the fact that V(:,1) is even
%and V(:,2) an odd function of x, i.e., f(x) = -f(-x).
%Thus, their product is an odd fucntion of x, the integral of an odd
%function from-A to +A is zero.
Vprod = conj(Eigen_Fun(:,2)).*(Eigen_Fun(:,1))*dy;
Vorth = sum(Vprod); %very small but not zero!
figure;
plot(y,Vnormalised1,'b',y,Vnormalised2,'r.',y,Vprod,'k')
legend('Ground state','First excited state','Orthogonality')
xlabel('$y$','Interpreter','latex')
ylabel('$\psi_i(y), i=1,2$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
%-------------

                            %%Eigenvalues%%%

%looking at the eignevalues
Eigen_Val = eig(H);
%plotting the eigenvalues
figure;
% %plot the entire spectrum of the eigenvalues
% plot(Eigen_Val)
% xlabel('Number of eigenvalues')
% ylabel('Eigenvalues')
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'FontSize',16)
%----
%plot the first four eigenvalues
%1st
h1 = plot(Eigen_Val(1),'ko','MarkerSize',8,'LineWidth',2);
set(h1, 'MarkerFaceColor', get(h1,'Color'));
% text(1.03,2.6,'$E_0$','Interpreter','latex','FontSize',24,'Color','black')
hold on
%2nd
h2 = plot(Eigen_Val(2),'bo','MarkerSize',8,'LineWidth',2);
set(h2, 'MarkerFaceColor', get(h2,'Color'));
% text(0.9,2.7,'$E_1$','Interpreter','latex','FontSize',24,'Color','blue')
%3rd
h3 = plot(Eigen_Val(3),'go','MarkerSize',8,'LineWidth',2);
set(h3, 'MarkerFaceColor', get(h3,'Color'));
% text(1.03,4.3,'$E_2$','Interpreter','latex','FontSize',24,'Color','green')
%4th
h4 = plot(Eigen_Val(4),'ro','MarkerSize',8,'LineWidth',2);
set(h4, 'MarkerFaceColor', get(h4,'Color'));
% text(1.03,4.6,'$E_3$','Interpreter','latex','FontSize',24,'Color','red')

set(gca,'xtick',[])
ylabel('Eigenvalues')
% ylim([1 5])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
%-------------

%Calculating the tunneling rate using two localised left and right
%the left localised numerical state/mode
phi_L = Vnormalised1+Vnormalised2/sqrt(2);
%the right localised numerical state/mode
phi_R = Vnormalised1-Vnormalised2/sqrt(2);

%Gaussian functions
DW_min = sqrt(b/(2*gamma));
%left Gaussian function
G_L = (1/pi)^(1/4).*exp(-(y+DW_min).^2);
%left Gaussian function
G_R = (1/pi)^(1/4).*exp(-(y-DW_min).^2);

%plotting
figure;
plot(y,-phi_L,'r',y,-phi_R,'b')
hold on
plot(y,G_L,'r.',y,G_R,'b.')
xlabel('y','Interpreter','latex')
ylabel('Ground State')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
legend('Left Localised Mode','Right Localised Mode',...
    'Left Gaussian Fucntion','Right Gaussian Fucntion')

%-------------




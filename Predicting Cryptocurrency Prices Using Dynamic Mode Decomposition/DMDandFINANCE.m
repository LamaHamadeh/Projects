
%This file uses Dynamic Mode Decomposition to predict the future prices of
%cryptocurrencies

clear all;
close all;

%% Define all interested cryptocurrencies
% Tickers = {'BTC-GBP' 'ETH-GBP' 'BNB-GBP' 'USDT-GBP' 'SOL-GBP' 'DOGE-GBP'...
%     'ADA-GBP'  'HEX-GBP'  'XRP-GBP'   'USDC-GBP' 'LUNA1-GBP' 'DOT-GBP'  ...
%     'AVAX-GBP' 'SHIB-GBP' 'MATIC-GBP' 'CRO-GBP'  'UNI1-GBP'  'LINK-GBP' ...
%     'LTC-GBP'  'ALGO-GBP' 'DAI-GBP'   'BCH-GBP'  'TRX-GBP'   'XLM-GBP'  ...
%     'MANA-GBP' 'ATOM-GBP' 'AXS-GBP'   'VET-GBP'  'SAND-GBP'  'FTT-GBP'  ...
%     'HBAR-GBP' 'FTM-GBP'  'FIL-GBP'   'ENJ-GBP'  'THETA-GBP' 'EGLD-GBP' ...
%     'ETC-GBP'  'XTZ-GBP'  'HNT-GBP'   'XMR-GBP'  'MIOTA-GBP' 'AAVE-GBP' ...
%     'GRT1-GBP' 'EOS-GBP'  'STX-GBP'   'CAKE-GBP' 'ONE1-GBP'  'LRC-GBP'  ...
%     'FLOW-GBP' 'BTT-GBP'}; %50

Tickers = {'ADA-GBP'  'HEX-GBP'  'XRP-GBP'   'USDC-GBP' 'LUNA1-GBP' 'DOT-GBP'  ...
    'AVAX-GBP' 'SHIB-GBP' 'MATIC-GBP' 'CRO-GBP'  'UNI1-GBP'  'LINK-GBP'...
    'LTC-GBP'  'ALGO-GBP' 'DAI-GBP'   'BCH-GBP'  'TRX-GBP'   'XLM-GBP'}; 

%retrieve historic financial data (prices) from yahoo! finance
% stock = hist_stock_data('01012019', '01012021', Tickers)
stock = hist_stock_data(now-10, now, Tickers)

%visualise the dynamics of crypto behaviour
% figure;
% % Bitcoin
% subplot(1,3,1)
% plot(stock(:,1).Open)
% xlabel('Days')
% ylabel('£')
% title('Bitcoin')
% axis square
% axis tight
% % Ethereum
% subplot(1,3,2)
% plot(stock(:,2).Open)
% xlabel('Days')
% ylabel('£')
% title('Ethereum')
% axis square
% axis tight
% % Dogecoin
% subplot(1,3,3)
% plot(stock(:,3).Open)
% xlabel('Days')
% ylabel('£')
% title('Dogecoin')
% axis square
% axis tight

%construct the dataset of all the historic prices
X = []; %NxM: N is the number of cryptocurrencies and M is the number of snapshots
for i = 1:length(Tickers)
    X = [X stock(:,i).Open];
end

%% Apply Dynamic Mode Decomposition (DMD) on prices data

%reconstruct it in a way that currencies are on the rows and time span is
%on the columns. The columns represent the snapshots of the prices.
X = X'; %N >> M
%----------------------------------------

%show the time dynamics of one of the cryptos
figure;
subplot(1,2,1)
plot(X(1,:))
xlabel('Time')
ylabel('Price £')
title('Before smoothing')
axis square
% axis([0 max(t) 1 1.15])
hold on
%----------------------------------------

%smooth the data (moving averages)
%calculating moving average on all data points
k = 3;% number of points within the sliding window of averaging
for i = 1:size(X,1)
    X(i,:) = movmean(X(i,:),k);
end
%show the smoothed time dynamics
subplot(1,2,2)
plot(X(1,:))
xlabel('Time')
ylabel('Price £')
title('After smoothing')
axis square 
% axis([0 max(t) 1 1.15])
%----------------------------------------

%define time variable (time window) and time step
dt = 1;
t = linspace(1,11,length(stock(:,1).Open));
%----------------------------------------

% create DMD matrices
X1 = X(:,1:end-1);
X2 = X(:,2:end);
%----------------------------------------

%SVD and rank truncation
[U,S,V] = svd(X1,'econ');
%check the modes that holds the most variance percentage
sigma = diag(S);
figure;
subplot(1,2,1)
plot(sigma/sum(diag(S)),'ko','LineWidth',2)
axis square
xlabel('$r$','interpreter','latex')
ylabel('$\displaystyle{\sigma_i/\sum_{i=1}^{N} \sigma_i}$','interpreter','latex')
title('Before Truncation')
hold on

%Define truncation rank
r = 4;
Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);
%show the percentage of information captured in each mode from the SVD
subplot(1,2,2)
sigmar = diag(Sr);
% semilogy(sigmar/trace(Sr),'ko','LineWidth',2)
% semilogy(sigmar/sum(diag(Sr)),'ko','LineWidth',2)
plot(sigmar/sum(diag(Sr)),'ko','LineWidth',2)
axis square
xlabel('$r$','interpreter','latex')
ylabel('$\displaystyle{\sigma_i/\sum_{i=1}^{r} \sigma_i}$','interpreter','latex')
title('After Truncation')
%----------------------------------------

%Build Atilde and DMD modes: least square matrix/linear operator
Atilde = Ur'*X2*Vr/Sr;
%eigen decomposition
[W,D] = eig(Atilde); %eigenvectors and eigenvalues
%DMD spatial modes 
Phi = X2*Vr/Sr*W; 
%visualise the DMD modes 
figure;
for i = 1:size(Phi,2)
    subplot(2,2,i)
    bar(real(Phi(:,i)))
end
%DMD eigenvalues spectra
lambda = diag(D);
omega = log(lambda)/dt; %frequency
%check the eignevalues that all of them imaginary
figure;
for ii = 1:length(omega(:,1))
    plot(real(omega(ii,:)),imag(omega(ii,:)),'r*','Linewidth',2)
    hold on
    xlabel('Real')
    ylabel('Imaginary')
end
%----------------------------------------

% reconstruction of modes/function in time
% compute DMD solution
x1 = X(:,1); %t=0
b = Phi\x1; %pseudo-inverse initial conditions
time_dynamics = zeros(r,length(t));
for iter = 1:length(t)
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
X_DMD = Phi * time_dynamics; %the full solution of the linear model generated by DMD

%show the behaviour of the DMD mode
figure;
%show the dynamics of the first crypto
subplot(1,2,1)
plot(t,X(1,:))
axis square
xlabel('Time')
ylabel('Price')
title('Actual Dynamics')
hold on
%show the dynamics of the reconstructrd DMD mode 
subplot(1,2,2)
plot(t,real(X_DMD(1,:)))
axis square
xlabel('Time')
ylabel('Price')
title('DMD Dynamics')

%plot time dynamics
figure;
for j = 1:length(time_dynamics(:,1))
    hold on
    subplot(2,2,j)
    plot(t,real(time_dynamics(j,:)),'r','LineWidth',2)
    xlabel('$t$','interpreter','latex')
    set(gca,'fontsize',24)
    axis square
    set(gca,'TickLabelInterpreter','latex')
    hold on
end
%----------------------------------------

% error estimation between DMD modes and the original data
for j=1:length(t)
    %average error
    E_avg(j)= sum(abs(X_DMD(:,j)-X(:,j)))/numel(X(:,j));
    %Root-mean-square error (RMSE): L2 error
    E_RMSE(j) = sqrt((sum(abs(X_DMD(:,j)-X(:,j)).^2))/(numel(X(:,j))));
    %maximum error
    E_max(j) = max(abs(X_DMD(:,j)-X(:,j)));
    %standard deviation
    stand(j) = std(X_DMD(:,j)-X(:,j));
end
figure;
%Average Error
subplot(2,2,1)
semilogy(t,E_avg,'k','LineWidth',2)
xlim([0 max(t)])
xlabel('$t$','interpreter','latex')
ylabel('$\displaystyle{\frac{\sum_{i=1}^{n}\left|X_{\mathrm{DMD}}-X\right|}{n}}$','interpreter','latex')
title('Average Error')
set(gca,'fontsize',14)
set(gca,'TickLabelInterpreter','latex')
hold on
%Root-mean-square error
subplot(2,2,2)
semilogy(t,E_RMSE,'k','LineWidth',2)
xlim([0 max(t)])
xlabel('$t$','interpreter','latex')
ylabel('$\displaystyle{\sqrt{\frac{\sum_{i=1}^{n}\left|X_{\mathrm{DMD}}-X\right|^2}{n}}}$','interpreter','latex')
title('Root-mean-square error')
set(gca,'fontsize',14)
set(gca,'TickLabelInterpreter','latex')
hold on
%Maximum Error
subplot(2,2,3)
semilogy(t,E_max,'k','LineWidth',2)
xlim([0 max(t)])
xlabel('$t$','interpreter','latex')
ylabel('$\displaystyle{max\left|X_{\mathrm{DMD}}-X\right|}$','interpreter','latex')
title('Maximum Error')
set(gca,'fontsize',14)
set(gca,'TickLabelInterpreter','latex')
hold on
%Standard deviation
subplot(2,2,4)
semilogy(t,stand,'k','LineWidth',2)
xlim([0 max(t)])
xlabel('$t$','interpreter','latex')
ylabel('$\displaystyle{std\left|X_{\mathrm{DMD}}-X\right|}$','interpreter','latex')
title('Standard deviation')
set(gca,'fontsize',14)
set(gca,'TickLabelInterpreter','latex')
hold on
%----------------------------------------

% short time prediction
t2 = linspace(0,20,20);
time_dynamics2 = zeros(r,length(t2));
for iter = 1:length(t2)
    time_dynamics2(:,iter) = (b.*exp(omega*t2(iter)));
end
X_DMD2 = Phi * time_dynamics2; %the full future solution of the linear model generated by DMD
%show the prediction
figure;
%original dynamics
subplot(1,3,1)
plot(t,X(1,:))
xlabel('$t$','interpreter','latex') 
ylabel('$X$','interpreter','latex')
axis square
set(gca,'fontsize',14)
set(gca,'TickLabelInterpreter','latex')
hold on
%DMD dynamics
subplot(1,3,2)
plot(t,real(X_DMD(1,:)))
xlabel('$t$','interpreter','latex')
ylabel('$X_\mathrm{DMD}$','interpreter','latex')
axis square
set(gca,'fontsize',14)
set(gca,'TickLabelInterpreter','latex')
hold on
%short time future prediction
subplot(1,3,3)
plot(t2,real(X_DMD2(1,:)))
xlabel('$t$','interpreter','latex')
ylabel('$X_\mathrm{DMD}$','interpreter','latex')
axis square
set(gca,'fontsize',14)
set(gca,'TickLabelInterpreter','latex')
%----------------------------------------

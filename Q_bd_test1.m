
clear
close all
format long

%% user input
epsr = 4; %material permittiivity
L = 1; %flim length
n = sqrt(epsr);
m_max = 12;
m_min = 0;
m0_max = 13;
m0_min = 0;
m = m_min:m_max;
m_bd = m_min:m_max;
m0_bd = m0_min:m0_max;

%% res freq

%unstructured
wr = m.*pi / n / L;
r = (n - 1) / (n + 1);
wi_FP = log(r) / n / L;
wi = wi_FP * ones(size(m));

%global bd
wr0_bd = m0_bd * pi / n / L;
wi0_bd = [-0.5493 -0.549 -0.477 -0.296 -0.176 -0.1 -0.0563 -0.0311 -0.0168 ...
    -0.00892 -0.00470 -0.00245 -0.00127 -0.00065]; %benchmark bound

%max violated local constraint (method 1) bd
wr_bd = m_bd * pi / n / L;
wi_bd = [-0.5493 -0.548 -0.477 -0.296 -0.171 -0.0980 -0.05250 -0.0275 ...
    -0.0135 -0.007 -0.0032 -0.0013 -0.0005250];
% wi_bd = [-0.5493 -0.545 -0.476 -0.296 -0.167 -0.0980 -0.0555 -0.0300 -0.0164];

%local optima
wr_loc = [1.57 3.6 5.88 8.34 10.57 13.14 15.33 17.96 20 22.73];
wi_loc = -[0.5494 0.5494 0.4517 0.2888 0.1912 0.1104 0.0680 0.0368 ...
    0.0212 0.0113];

%analytical bound
wr_bd_ana = 0:0.01:1;
wi_bd_ana = -3 * exp(-n^2*L^2*wr_bd_ana.^2) ./ (exp(1) * L) ...
    ./ (1 + 3*n^2 + n^2*L^2*wr_bd_ana.^2);

t = [4 4 25 70 160]; %computation time (seconds)
DOF = [50 50 100 100 125]; %degrees of freedom
ND_max = [10 10 15 30 35]; %max ND needed for reasonable convergence 

%% Q factor

Q = wr ./ (-2*wi);
Q_loc = wr_loc ./ (-2*wi_loc);
Q_bd0 = wr0_bd ./ (-2*wi0_bd);
Q_bd = wr_bd ./ (-2*wi_bd);
Q_bd_ana = wr_bd_ana ./ (-2*wi_bd_ana); 

%% plot

set(0,'DefaultLineLineWidth', 2)
set(0,'defaultAxesFontSize', 16)
set(0, 'DefaultAxesBox', 'on')
figure
plot(wr, Q, '-ko')
hold on
% plot(wr_loc, Q_loc, '-bo')
% plot(wr0_bd, Q_bd0, '-ro')
% plot(wr_bd, Q_bd, '--go')
plot(wr_bd_ana,Q_bd_ana,'m')
% plot(wr_bd, Q_bd_f, '-go')
xlabel('Frequency, \omega')
ylabel('Quality factor')
title('FP Cavity Quality Factor, n = 2')
axis tight
xlim([0 10])
ylim([0 60])
% ylim([0,1100])
% legend({'Unstructured','local optima','global bound','most violated constraint','analytical bound'},'location','best')
legend({'Unstructured','analytical bound'})

figure(2)
semilogy(wr, Q, '-k');
hold on
semilogy(wr_bd_ana,Q_bd_ana,'m')
% semilogy(wr_loc,Q_loc,'-bo')
semilogy(wr0_bd,Q_bd0,'-ro')
semilogy(wr_bd,Q_bd,'--go')
% semilogy(wr_bd_ana,Q_bd_ana,'m')
% semilogy(wr_bd,Q_bd_f,'-go')
xlabel('Frequency, \omega')
ylabel('Quality factor')
title('FP Cavity Quality Factor, n = 2')
axis tight
% ylim([0,1100])
legend({'Unstructured','analytical bound','global bound','most violated constraint'},'location','best')

% figure(3)
% semilogy(wr, Q, '-ko')
% hold on
% semilogy(wr_bd_ana,Q_bd_ana,'m')
% xlabel('Frequency, \omega')
% ylabel('Quality factor')
% title('FP Cavity Quality Factor, n = 2')
% legend({'Unstructured','Analytical Bound'})
% axis tight
% xlim([0 3*pi])
% ylim([0 60])
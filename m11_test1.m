clear
close all
format long
addpath(genpath('function'))

%% user input

tic
epsr = 4; %material permittiivity
L = 1; %flim length
L_perb = 0.0; %perturbation (if fourier space doesnt find feasible solution)
L = L + L_perb;
Nx = 205; %resolution
Nf = 150;
ND = 205; %num local power constraints
m = 11;
z = linspace(0,L,Nx);

%% resonant frequencies

n = sqrt(epsr); %refractive index
wr = m*pi/n/L;
wi = -0.0013;%<---- MAX IMAG(W)
w = wr + 1j*wi;
r = (n-1) / (n+1);
wi_FP = log(r) / n / L;
wavelength = 2*pi / (wr);
fprintf('wr = %4.5f\n', wr)
fprintf('wi guess = %4.5f\n', wi)
fprintf('wavelength = %4.2f\n', wavelength)

%% init

xi = get_xi(epsr);
G0 = cal_G0_1D(z, w);
Ngrid = length(G0);
xi_Mat = xi*eye(Ngrid);
S = G0 + xi_Mat; % constraint
% U = x2f(z,Nf,L*(1+L_perb)); % fourier space matrices
U = eye(Nx,Nx); % real space
% [nrow,ncol] = size(U);
A = - (imag(w'*G0) + imag(xi_Mat)); %objective function

%% D matrices

D = cell(1,ND);
D{1} = (eye(size(S)) * 0);
D{2} = (eye(size(S)) * 1);
fmax = zeros(1,ND);
popt = cell(1,ND);
ortho = cell(1,ND);

%% ImG0 bound

[fmax(1),popt{1}] = bound_Ds(S,{D{1}},A,epsr,U); % ImG0 bound
[fmax(2),popt{2}] = bound_Ds(S,D,A,epsr,U); % ReG0 bound

%% iterate D matrix bounds

for i = 3:ND
    einc = zeros(size(G0,1),1);
    D{i} = get_Dopt(popt{i-1},S,einc);
    [fmax(i),popt{i}] = bound_Ds(S,D,A,epsr,U);
    fprintf(' D matrix: %d / %d, fmax = %s \n', i, ND, fmax(i))
end

t = toc;
fprintf('Elapsed Time = %4.2f seconds \n',t)

%% result

wi_bound = fmax;
set(0,'DefaultLineLineWidth', 2)
set(0,'defaultAxesFontSize', 16)
set(0, 'DefaultAxesBox', 'on')
figure
x = 1:ND;
hold on
plot(x,wi_bound,'-x')
plot(x,wi*ones(1,ND))
plot(x,wi_FP*ones(1,ND))
xlim([1,ND])
xlabel('Number of D matrices, ND')
ylabel('Im\omega')
legend({'upper bound','initial guess','FP cavity'})
title(tName({'m','\epsilon','L','Nx','Im{\omega}','\delta t (sec)'},{m,epsr,L,Nx,wi,t},''))


% figure
% for i = 3:ND
%     plot(real(diag(D{i})))
%     hold on
% end

check = ones(1,ND);
for i = 3:ND
    if fmax(i) >= fmax(i-1)
        check(i) = 0;
    end
end
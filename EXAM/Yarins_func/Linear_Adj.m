clc;
clear all;
close all;


%--------------------------------------------------------------------------
%   Observations and redundancy
%--------------------------------------------------------------------------

x = [17.95 14.99 13.91 17.01 15.94 17.96 9.97 19.04 14.02]';
y =[56.94 58.98 50.93 50.01 57.03 52.96 53.01 50.08 57.97]';
%Vector of observations (FILL!!!)
L = [12.82 13.12 7.43 7.83 12.13 10.15 7.48 8.45 9.47]';

%Number of observations
no_n = length(L);

%Number of unknowns (FILL!!!)
no_u = 6;

%Redundancy
r = no_n - no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations (FILL!!!)
%s_L = [0.02 0.02 0.02 0.02 0.02];
%S_LL = diag(s_L.^2);

%Theoretical reference standard deviation
sigma_0 = 1;        %a priori

%Cofactor matrix of the observations
%Q_LL = (1 / sigma_0^2) * S_LL;
Q_LL = eye(9);

%Weight matrix
%P = inv(Q_LL);
P = eye(9);
%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Design matrix
%x = [1 2 3 4 5]';
A = [ones(9,1) x y x.*y x.^2 y.^2]

%Normal matrix
N = A'*P*A;

%Vector of right hand side of normal equations
n = A'*P*L;
    
%Inversion of normal matrix / Cofactor matrix of the unknowns
Q_XX = inv(N);
    
%Solution of normal equation
X_hat = Q_XX*n;
    
%Estimated unknown parameters
a = X_hat(1)
b = X_hat(2)
c = X_hat(3)

%Vector of residuals
v = A * X_hat - L;

%Objective function
vTPv = v' * P * v;

%Vector of adjusted observations
L_hat = L + v;

%Final check
check = A*X_hat-L_hat

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);    %a posteriori
r
%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_XX;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_XX*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));

%Cofactor matrix of the residuals
Q_vv = Q_LL - Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));



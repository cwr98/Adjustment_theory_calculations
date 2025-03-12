clc;
clear all;
close all;


%--------------------------------------------------------------------------
%   Observations and redundancy
%--------------------------------------------------------------------------

%Vector of observations (FILL!!!)
L = [10.509 5.360 -8.523 -7.348 -3.167 15.881]';
cons = [100 0 0 -100 0 100]';
L_mod = L+cons;

%Number of observations
no_n = length(L);

%Number of unknowns (FILL!!!)
no_u = 3;

%Redundancy
r = no_n - no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations (FILL!!!)
s_L = [6 4 5 3 4 12]/1000;
S_LL = diag(s_L.^2);

%Theoretical reference standard deviation
sigma_0 = 1;        %a priori

%Cofactor matrix of the observations
Q_LL = (1 / sigma_0^2) * S_LL;
%Q_LL = eye(no_n);

%Weight matrix
P = inv(Q_LL);
%P = eye(no_n);
%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Design matrix
%x = [1 2 3 4 5]';
A = [1 0 0;
    -1 1 0;
    0 -1 1;
    0 0 -1;
    -1 0 1;
    0 1 0];

%Normal matrix
N = A'*P*A;

%Vector of right hand side of normal equations
n = A'*P*L_mod;
    
%Inversion of normal matrix / Cofactor matrix of the unknowns
Q_XX = inv(N);
    
%Solution of normal equation
X_hat = Q_XX*n;
    
%Estimated unknown parameters
a = X_hat(1)
b = X_hat(2)
c = X_hat(3)

%Vector of residuals
v = A * X_hat - L_mod;

%Objective function
vTPv = v' * P * v;

%Vector of adjusted observations
L_hat = L + v;

%Final check
check = A*X_hat-L_hat-cons

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);    %a posteriori

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



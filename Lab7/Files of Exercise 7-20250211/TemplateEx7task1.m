%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%    Exercise 7: Adjustment Calculation - part II  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 09, 2018
%   Last changes   : December 14, 2022
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Task 1: Adjustment of a straight line
%--------------------------------------------------------------------------
x=[1.0; 2.0; 3.0; 4.0];
y=[0.1; 1.1; 1.8; 2.4];

%--------------------------------------------------------------------------
%   Observations and redundancy
%--------------------------------------------------------------------------
%Vector of observations 
L = y;

%Number of observations
no_n = length(L);  

%Number of unknowns
no_u = 2;

%Redundancy
r = no_n-no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_L=[2 1 4 2];

S_LL = diag(s_L.^2);

%Theoretical reference standard deviation
sigma_0 = 1;      %a priori

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Design matrix
A = [x ones(no_n,1)];
    
%Normal matrix
N = A'*P*A;
        
%Vector of right hand side of normal equations
n = A'*P*L;
    
%Inversion of normal matrix / Cofactor matrix of the unknowns
Q_XX = inv(N);
    
%Solution of normal equation
X_hat = Q_XX*n;
    
%Estimated unknown parameters
a_hat=X_hat(1);
b_hat=X_hat(2);
 
%Vector of residuals
v = A*X_hat-L;

%Objective function
vTPv = v'*P*v;

%Vector of adjusted observations
L_hat = L+v;

%Final check
max(abs(L_hat-(a_hat*x+b_hat)))
  
%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);       %a posteriori

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
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));

disp('adjusted unknowns and std')
[X_hat s_X]

disp('residuals and adj. obs. and their std.')
[v s_v L_hat s_L_hat]



plot(x,y,'.','markersize',15)
hold on
plot(x,a_hat*x+b_hat)



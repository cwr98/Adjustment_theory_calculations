clc;
clear all;
close all;


%--------------------------------------------------------------------------
%   Observations and redundancy
%--------------------------------------------------------------------------

%Vector of observations (FILL!!!)
L = [3.0 1.5 0.2]';

%Number of observations
no_n = length(L);

%Number of unknowns (FILL!!!)
no_u = 2;

%Number of constraints
no_b = 1;

%Redundancy
r = no_n-no_u+no_b;  

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations (FILL!!!)
s_L = [4 2 1]./100;
S_LL = diag(s_L.^2);

%Theoretical reference standard deviation
sigma_0 = 1;        %a priori

%Cofactor matrix of the observations
Q_LL = (1 / sigma_0^2) * S_LL;
%Q_LL = eye(7);

%Weight matrix
P = inv(Q_LL);
%P = eye(7);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Design matrix
%x = [1 2 3 4 5]';
A = [1 1;
    2 -1;
    1 -1];

 %Design matrix C with the elements from the Jacobian matrix J
C = [0.9 -1]; %for the first constraint

 %Normal matrix
 N = A'*P*A;
 
 %Extended normal matrix
 N_ext = [N C';
     C 0];
 
 %Vector of right hand side of normal equations
 n = A'*P*L;

 % CX=c
 c = 0.0;
 
 %Extended vector of right hand side of normal equations
 n_ext = [n; c];

 %Inversion of normal matrix / Cofactor matrix of the unknowns
 Q = inv(N_ext);
 Q_XX = Q(1:no_u,1:no_u);

 %Solution of the normal equations
 x_solution = Q*n_ext;
    
%Estimated unknown parameters
x = x_solution(1)
y = x_solution(2)
k = x_solution(3)

X_hat = x_solution(1:no_u);

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



%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%    Exercise 8: Adjustment Calculation - part III  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 09, 2018
%   Last changes   : January 04, 2023
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Task 1 - Non-linear equation system
%--------------------------------------------------------------------------
disp('Task 1 - Non-linear adjustment problem!')

%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Vector of observations
L = [-4; 8; 7.7;];

%Number of observations
no_n = length(L);

%Initial values for the unknowns
x = 2;
y = 2;

%Vector of initial values for the unknowns
X_0 = [x y]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r =  no_n - no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
%S_LL = not given!

%Theoretical standard deviation
%sigma_0 =      %a priori

%Cofactor matrix of the observations
Q_LL = eye(no_n);

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-12;
delta = 10^-12;
max_x_hat = Inf;
Check2 = Inf;

%Number of iterations
iteration = 0;

while max_x_hat>epsilon || Check2>delta
    
     %Observations as functions of the approximations for the unknowns
    L_0 = [x+y-2*y^2; x^2+y^2; 3*x^2-y^2];

     
     %Vector of reduced observations
     l = L - L_0;
    
     %Design matrix with the elements from the Jacobian matrix J
     A = [1 1-4*y; 2*x 2*y; 6*x -2*y];
    
     %Normal matrix
     N = A'*P*A;
     
     %Vector of right hand side of normal equations
     n = A'*P*l;
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx = inv(N);
    
     %Solution of the normal equations
     x_hat = Q_xx*n;
       
     %Update
     X_0 =  X_0+x_hat;
    
     x = X_0(1);
     y = X_0(2);

    
     %Check 1
    max_x_hat = max(abs(x_hat));
     
     %Vector of residuals
     v = A*x_hat-l;
 
     %Vector of adjusted observations
     L_hat = L+v;
    
     %Objective function
     vTPv = v'*P*v;
    

     %Check 2
     Check2 = max(abs(L_hat-[x+y-2*y^2; x^2+y^2; 3*x^2-y^2]));
    %this check not so necessery, its suffieciet to just check dX.

     %Update number of iterations
     iteration = iteration+1;
  
end

%rest copy from previous exercise:

%Empirical reference standard deviation
%s_0 = 

%VC matrix of adjusted unknowns
%S_XX_hat = 

%Standard deviation of the adjusted unknowns
%s_X = 

%Cofactor matrix of adjusted observations
%Q_LL_hat = 

%VC matrix of adjusted observations
%S_LL_hat = 

%Standard deviation of the adjusted observations
%s_L_hat = 

%Cofactor matrix of the residuals
%Q_vv = 

%VC matrix of residuals
%S_vv = 

%Standard deviation of the residuals
%s_v = 



%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%    Exercise 12: Adjustment Calculation - part VII  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 12, 2018
%   Last changes   : January 31, 2023
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
format long g;

%--------------------------------------------------------------------------
%   Task 1
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Load files
% dir =                %[gon]
% coord =              %[m]     %Error-free



%Vector of observations
% L = 

%Number of observations
no_n = length(L);

%Initial values for the unknowns
% x3 = 
% y3 = 
% w3 = 

%Vector of initial values for the unknowns
% X_0 = 

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n-no_u;  

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
% s_dir = 

%VC Matrix of the observations
% S_LL =   

%Theoretical standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-5;
delta = 10^-13;
max_x_hat = Inf;
Check2 = Inf;

%Number of iterations
iteration = 0;

while max_x_hat>epsilon || Check2>delta            
    
     %Observations as functions of the approximations for the unknowns

     %Vector of reduced observations
%      l = 
    
     %Design matrix with the elements from the Jacobian matrix J

    
     %Normal matrix
%      N = 
     
     %Vector of right hand side of normal equations
%      n = 
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
%      Q_xx = 
    
     %Solution of the normal equations
%      x_hat = 
       
     %Update
%      X_0 =  

    
     %Check 1
%      max_x_hat =
     
     %Vector of residuals
%      v = 
 
     %Vector of adjusted observations
%      L_hat =
    
     %Function
%      vTPv = 
    
     %Functional relationships without the observations
%      phi_X_hat = 

     %Check 2
%      Check2 = 
    
     %Update number of iterations
     iteration = iteration+1;
  
end

if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end





%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));        %[m]

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

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

    


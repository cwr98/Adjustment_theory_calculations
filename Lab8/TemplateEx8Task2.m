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
%   Task 2 - Copper 
%--------------------------------------------------------------------------
disp('Task 2 - Non-linear adjustment problem!')

%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------


%Vector of observations
%L = 

%Number of observations
%no_n = 

%Initial values for the unknowns


%Number of unknowns
%no_u = 

%Redundancy
%r =   

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------


%VC Matrix of the observations
%S_LL = 

%Theoretical standard deviation
%sigma_0 =     %a priori

%Cofactor matrix of the observations
%Q_LL = 

%Weight matrix
%P = 

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
%epsilon = 
%delta = 
%max_x_hat = 

%Number of iterations
iteration = 0;

%while      
    
     %Observations as functions of the approximations for the unknowns
     %L_0
	 
     
     %Vector of reduced observations
     %l = 
    
     %Design matrix with the elements from the Jacobian matrix J
     %A = 
    
     %Normal matrix
     %N = 
     
     %Vector of right hand side of normal equations
     %n = 
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
     %Q_xx = 
    
     %Solution of the normal equations
     %x_hat = 
       
     %Update
     
    
     %Check 1
     %max_x_hat = 
     
     %Vector of residuals
     %v = 
 
     %Vector of adjusted observations
     %L_hat = 
    
     %Objective function
     %vTPv = 
    


     %Check 2
     %Check2 = 
    
     %Update number of iterations
     iteration = iteration+1;
  
%end



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


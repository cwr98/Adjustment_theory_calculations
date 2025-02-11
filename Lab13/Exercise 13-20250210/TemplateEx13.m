%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%   Exercise 13: Adjustment Calculation - part VIII  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 15, 2018
%   Last changes   : February 08, 2023
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
% Observations and initial values for unknowns
%--------------------------------------------------------------------------
% Load all files
%dist = 
%dir = 
%fixedpoint = 
%newpoint = 

% Vector of observations
%L = 

% Gauss-Krueger coordinates for control points [m]
for i=1:size(fixedpoint,1)
    eval(['y' num2str(fixedpoint(i,1)) '=' num2str(fixedpoint(i,2)) ';']);
    eval(['x' num2str(fixedpoint(i,1)) '=' num2str(fixedpoint(i,3)) ';']);
end

% Gauss-Krueger coordinates for new points [m]
for i=1:size(newpoint,1)
    eval(['y' num2str(newpoint(i,1)) '=' num2str(newpoint(i,2)) ';']);
    eval(['x' num2str(newpoint(i,1)) '=' num2str(newpoint(i,3)) ';']);
end

% Initial values for orientation unknowns
%w1 = 
%w6 = 
%w9 = 
%w15 = 

% Initial values for unknowns
%X_0 = 

% Number of observations
no_n = length(L);

% Number of unknowns
no_u = length(X_0);
 
% Redundancy
r = no_n-no_u;


%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
% VC Matrix of the observations

%S_LL = 

% Theoretical standard deviation
sigma_0 = 1;

% Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

% Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
% break-off condition
epsilon = 10^-5;
delta = 10^-15;
max_x_hat = Inf;
Check2=Inf;

% Number of iterations
iteration = 0;

% Initialising A
A = zeros(no_n,no_u);

% Iteration
while max_x_hat>epsilon || Check2>delta

	% Distances
	%L_0 = 
	
	% Directions
	%L_0=

    
    % Vector of reduced observations
    l = L-L_0';

    % Design matrix with the elements from the Jacobian matrix J
	%A =
  
    % Normal matrix
    N = A'*P*A;

    % Vector of right hand side of normal equations
    n = A'*P*l;

    % Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx = inv(N);
    
    % Solution of normal equations
    x_hat = Q_xx*n;
    
    % Adjusted unknowns
    X_hat = X_0+x_hat;
    
    % Update
   
    
    % Check 1
    max_x_hat = max(abs(x_hat));
    
    % Vector of residuals
    v = A*x_hat-l;
 
    % Vector of adjusted observations
    L_hat = L+v;
    
    % Function
    vTPv = v'*P*v;
    
    % Functional relationships 
    %phi_X_hat = 
	 
    % Check 2
    %Check2 = 
    
    % Update number of iterations
    iteration = iteration+1;

end

% Convert to [gon]





% Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

% VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

% Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));


% Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

% VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

% Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));


% Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

% VC matrix of residuals
S_vv = s_0^2*Q_vv;              

% Standard deviation of the residuals
s_v = sqrt(diag(S_vv));





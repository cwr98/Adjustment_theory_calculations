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
clearvars;
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
 dir = load('directions.txt');               %[gon]
 coord = load('Points.txt');             %[m]     %Error-free
 for i = 1:size(coord,1)
    eval(['y' num2str(coord(i,1)) '=' num2str(coord(i,2)) ';']);
    eval(['x' num2str(coord(i,1)) '=' num2str(coord(i,3)) ';']);
 end

%Vector of observations
L = dir(:,3)*pi/200; %[rad]

%Number of observations
no_n = length(L);

%Initial values for the unknowns
x3 = 250;
y3 = 500;
w3 = 1;

%Vector of initial values for the unknowns
X_0 = [x3 y3 w3]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n-no_u;  

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
s_dir = 0.001 * pi / 200; 

%VC Matrix of the observations
S_LL = s_dir^2*eye(no_n);  

%Theoretical standard deviation
sigma_0 = 0.001;     %a priori

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
     L_0(1,1) = direction(y3, x3, y1, x1, w3);
     L_0(2,1) = direction(y3, x3, y2, x2, w3);
     L_0(3,1) = direction(y3, x3, y4, x4, w3);
     L_0(4,1) = direction(y3, x3, y5, x5, w3);
     L_0(5,1) = direction(y3, x3, y6, x6, w3);
     %Vector of reduced observations
     l = L - L_0;
    
     %Design matrix with the elements from the Jacobian matrix J
     A(1,1) = dr_dx_from(y3,x3,y1,x1);
     A(1,2) = dr_dy_from(y3,x3,y1,x1);
     A(1,3) = -1;
     A(2,1) = dr_dx_from(y3,x3,y2,x2);
     A(2,2) = dr_dy_from(y3,x3,y2,x2);
     A(2,3) = -1;
     A(3,1) = dr_dx_from(y3,x3,y4,x4);
     A(3,2) = dr_dy_from(y3,x3,y4,x4);
     A(3,3) = -1;
     A(4,1) = dr_dx_from(y3,x3,y5,x5);
     A(4,2) = dr_dy_from(y3,x3,y5,x5);
     A(4,3) = -1;
     A(5,1) = dr_dx_from(y3,x3,y6,x6);
     A(5,2) = dr_dy_from(y3,x3,y6,x6);
     A(5,3) = -1;
    
     %Normal matrix
     N = A' * P * A;
     
     %Vector of right hand side of normal equations
    n = A' * P * l;
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx = inv(N);
    
     %Solution of the normal equations
     x_hat = Q_xx * n;
       
     %Update
     X_hat = X_0 + x_hat;
     X_0 = X_hat;
     x3 = X_hat(1);
     y3 = X_hat(2);
     w3 = X_hat(3); 

    
     %Check 1
     max_x_hat = max(abs(x_hat));
     
     %Vector of residuals
     v = A * x_hat - l;
 
     %Vector of adjusted observations
     L_hat = L + v;
    
     %Objective function
     vTPv = v' * P * v; 
    
     %Functional relationships without the observations
     phi_X_hat(1,1) = direction(y3, x3, y1, x1, w3);
     phi_X_hat(2,1) = direction(y3, x3, y2, x2, w3);
     phi_X_hat(3,1) = direction(y3, x3, y4, x4, w3);
     phi_X_hat(4,1) = direction(y3, x3, y5, x5, w3);
     phi_X_hat(5,1) = direction(y3, x3, y6, x6, w3);

     %Check 2
    Check2 = max(abs(L_hat - phi_X_hat));
    
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


%-------------------------------------------------------------------------
%   Task 2
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%-------------------------------------------------------------------------
%Load files
dir = load('directions.txt');               
%[gon]
coord = load('Points.txt');             
%[m]     %Error-free
for i = 1:size(coord,1)
eval(['y' num2str(coord(i,1)) '=' num2str(coord(i,2)) ';']);
eval(['x' num2str(coord(i,1)) '=' num2str(coord(i,3)) ';']);
end
%Vector of observations
L_ang(1) = (dir(3,3)-dir(2,3))*pi/200; %[rad]
L_ang(2) = (dir(4,3)-dir(3,3))*pi/200; %[rad]
L_ang(3) = (dir(5,3)-dir(4,3))*pi/200; %[rad]
L_ang(4) = (dir(1,3)-dir(5,3))*pi/200; %[rad]
L_ang = L_ang';
%Number of observations
no_n = length(L_ang);
%Initial values for the unknowns
x3 = 200;
y3 = 500;
w3 = 0;
%Vector of initial values for the unknowns
X_0_ang = [x3 y3]';
%Number of unknowns
no_u = length(X_0_ang);
%Redundancy
r = no_n-no_u;  

%-------------------------------------------------------------------------
%  Stochastic model
%-------------------------------------------------------------------------
s_dir = 0.001 * pi / 200;
    
%VC Matrix of the observations
S_LL_dir = s_dir^2*eye(5);
F = [0 -1 1 0 0;           
 0 0 -1 1 0;
 0 0 0 -1 1;
 1 0 0 0 -1];
S_LL = F * S_LL_dir * F';   
% Variance-covariance propagation
%Theoretical standard deviation
sigma_0 = 1e-5;     
%a priori
%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;
%Weight matrix
P = inv(Q_LL);
%-------------------------------------------------------------------------
%  Adjustment
%-------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-5;
delta = 10^-13;
max_x_hat = Inf;
Check2 = Inf;
%Number of iterations
iteration = 0;
while max_x_hat>epsilon || Check2>delta            
    %Observations as functions of the approximations for the unknowns
     L_0_ang(1,1) = direction(y3, x3, y4, x4, w3) - direction(y3, x3, y2, x2, w3);
     L_0_ang(2,1) = direction(y3, x3, y5, x5, w3) - direction(y3, x3, y4, x4, w3);
     L_0_ang(3,1) = direction(y3, x3, y6, x6, w3) - direction(y3, x3, y5, x5, w3);
     L_0_ang(4,1) = direction(y3, x3, y1, x1, w3) - direction(y3, x3, y6, x6, w3);
    %Vector of reduced observations
     l = L_ang - L_0_ang;
    %Design matrix with the elements from the Jacobian matrix J
     A_ang(1,1) = der_ang_x(y3,x3,y4,x4,y2,x2);
     A_ang(1,2) = der_ang_y(y3,x3,y4,x4,y2,x2);
     A_ang(2,1) = der_ang_x(y3,x3,y5,x5,y4,x4);
     A_ang(2,2) = der_ang_y(y3,x3,y5,x5,y4,x4);
     A_ang(3,1) = der_ang_x(y3,x3,y6,x6,y5,x5);
     A_ang(3,2) = der_ang_y(y3,x3,y6,x6,y5,x5);
     A_ang(4,1) = der_ang_x(y3,x3,y1,x1,y6,x6);
     A_ang(4,2) = der_ang_y(y3,x3,y1,x1,y6,x6);
    %Normal matrix
     N = A_ang' * P * A_ang;
    %Vector of right hand side of normal equations
     n = A_ang' * P * l;
    %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx = inv(N);
    %Solution of the normal equations
     x_hat = Q_xx * n;
    %Update
     X_hat_ang = X_0_ang + x_hat;
     X_0_ang = X_hat_ang;
     x3 = X_hat_ang(1);
     y3 = X_hat_ang(2);
    %Check 1
     max_x_hat = max(abs(x_hat));
    %Vector of residuals
     v_ang = A_ang * x_hat - l;
    %Vector of adjusted observations
     L_hat_ang = L_ang + v_ang;
    %Objective function
     vTPv = v_ang' * P * v_ang;
    %Functional relationships without the observations
     phi_X_hat_ang(1,1) = direction(y3, x3, y4, x4, w3) - direction(y3, x3, y2, x2, w3);
     phi_X_hat_ang(2,1) = direction(y3, x3, y5, x5, w3) - direction(y3, x3, y4, x4, w3); 
     phi_X_hat_ang(3,1) = direction(y3, x3, y6, x6, w3) - direction(y3, x3, y5, x5, w3);
     phi_X_hat_ang(4,1) = direction(y3, x3, y1, x1, w3) - direction(y3, x3, y6, x6, w3);
    %Check 2
     Check2 = max(abs(L_hat_ang - phi_X_hat_ang));
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
s_X_ang = sqrt(diag(S_XX_hat));        %[m]
%Cofactor matrix of adjusted observations
Q_LL_hat = A_ang*Q_xx*A_ang';
%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;
%Standard deviation of the adjusted observations
s_L_hat_ang = sqrt(diag(S_LL_hat));     
%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;
%VC matrix of residuals
S_vv = s_0^2*Q_vv;
%Standard deviation of the residuals
s_v_ang = sqrt(diag(S_vv));         
table(X_hat(1:2,:), X_hat_ang, s_X(1:2,:), s_X_ang, 'RowNames',{'x3' 'y3'}, 
'VariableNames',["X_hat", "X_hat_ang", "s_X", "s_X_ang"])

table(L, v, L_hat, s_v, s_L_hat, 'RowNames', {'r31' 'r32' 'r34' 'r35' 'r36'})
table(L_ang, v_ang, L_hat_ang, s_v_ang, s_L_hat_ang, 'RowNames',{'alpha' 'beta' 
'gamma' 'delta'})

clc;
clear all;
close all;
format long g;

x1 = 321.052;
y1 = 682.415;

x2 = 310.527;
y2 = 203.526;

x4 = 506.222;
y4 = 251.992;

x5 = 522.646;
y5 = 420.028;

x6 = 501.494;
y6 = 594.553;

x7 = 219.089;
y7 = 485.959;

%Vector of observations
L = [206.9094 46.5027 84.6449 115.5251 155.5891]'.*pi/200;

%Number of observations
no_n = length(L);

%Initial values for the unknowns
x3 = 300; 
y3 = 500;
w = 1;

%Vector of initial values for the unknowns
X_0 = [x3 y3 w]';

%Number of unknowns
no_u = length(X_0);

%Number of constraints
no_b = 1;

%Redundancy
r = no_n-no_u+no_b;  

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_r = (1/1000)*pi/200;
s_l = [s_r s_r s_r s_r s_r];
S_LL = diag(s_l.^2);  

%Theoretical standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_LL = (1/sigma_0^2)*S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-5;
delta = 10^-8;
max_x_hat = Inf;
Check2 = Inf;

%Number of iterations
iteration = 0;

while max_x_hat>epsilon || Check2>delta            
    
    %Observations as functions of the approximations for the unknowns
    L_0 = [direction_rad(x3,y3,x1,y1,w) direction_rad(x3,y3,x2,y2,w) direction_rad(x3,y3,x4,y4,w) direction_rad(x3,y3,x5,y5,w) direction_rad(x3,y3,x6,y6,w)]';
     
     %Vector of reduced observations
     l = L-L_0;
    
     %Design matrix A with the elements from the Jacobian matrix J
     A = [direction_der(x3,y3,x1,y1,w,'x1') direction_der(x3,y3,x1,y1,w,'y1') -1;
         direction_der(x3,y3,x2,y2,w,'x1') direction_der(x3,y3,x2,y2,w,'y1') -1;
         direction_der(x3,y3,x4,y4,w,'x1') direction_der(x3,y3,x4,y4,w,'y1') -1;
         direction_der(x3,y3,x5,y5,w,'x1') direction_der(x3,y3,x5,y5,w,'y1') -1;
         direction_der(x3,y3,x6,y6,w,'x1') direction_der(x3,y3,x6,y6,w,'y1') -1];

     
     %Design matrix C with the elements from the Jacobian matrix J
    C = [distance_der(x3,y3,x7,y7,'x1') distance_der(x3,y3,x7,y7,'y1') 0]; %for the first constraint
    
     %Normal matrix
     N = A'*P*A;
     
     %Extended normal matrix
     N_ext = [N C';
         C 0];
     
     %Vector of right hand side of normal equations
     n = A'*P*l;

     %vector of misclosure
     w = distance_(x3,y3,x7,y7)-25.0;
     
     %Extended vector of right hand side of normal equations
     n_ext = [n; -w];
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q = inv(N_ext);
     Q_xx = Q(1:no_u,1:no_u);
    
     %Solution of the normal equations
     x_solution = Q*n_ext;
     x_hat = x_solution(1:no_u);
       
     %Update
     X_hat = X_0+x_hat;
     X_0 = X_hat;

     x3=X_hat(1);
     y3=X_hat(2);
     w=X_hat(3);
    

     
     %Lagrange multiplier
     k = x_solution(end);
    
     %Check 1
     max_x_hat = max(abs(x_hat));
     
     %Vector of residuals
     v = A*x_hat-l;
 
     %Vector of adjusted observations
     L_hat = L+v;
    
     %Function
     vTPv = v'*P*v;
    
     %Functional relationships without the observations
     phi_X_hat = [direction_rad(x3,y3,x1,y1,w) direction_rad(x3,y3,x2,y2,w) direction_rad(x3,y3,x4,y4,w) direction_rad(x3,y3,x5,y5,w) direction_rad(x3,y3,x6,y6,w)]';

     %Check 2
     Check2 = max(abs(L_hat-phi_X_hat));
    
     %Update number of iterations
     iteration = iteration+1;
  
end

if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

x3
y3
w


%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));        

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

    


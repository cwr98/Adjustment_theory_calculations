
clc;
clear all;
close all;

%fixed values
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


%Vector of observations (distances)
L = [206.9094 46.5027 84.6449 115.5251 155.5891]'.*pi/200; %rad

%Number of observations
no_n = length(L);

%Initial values for the unknowns  
x3 = 200;
y3 = 460;
w = 1;

%Vector of initial values for the unknowns
X_0 = [x3 y3 w]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n - no_u; 

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_r = 1/1000; %gon
s_L = [s_r s_r s_r s_r s_r].*pi/200; %rad
S_LL = diag(s_L.^2);

%Theoretical standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_LL = (1/sigma_0^2)*S_LL;
%Q_LL =eye(9);

%Weight matrix
P = inv(Q_LL);
%P= eye(9);

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
    
     %Design matrix with the elements from the Jacobian matrix J
     A = [direction_der(x3,y3,x1,y1,w,'x1') direction_der(x3,y3,x1,y1,w,'y1') direction_der(x3,y3,x1,y1,w,'omega');
         direction_der(x3,y3,x2,y2,w,'x1') direction_der(x3,y3,x2,y2,w,'y1') direction_der(x3,y3,x2,y2,w,'omega');
         direction_der(x3,y3,x4,y4,w,'x1') direction_der(x3,y3,x4,y4,w,'y1') direction_der(x3,y3,x4,y4,w,'omega');
         direction_der(x3,y3,x5,y5,w,'x1') direction_der(x3,y3,x5,y5,w,'y1') direction_der(x3,y3,x5,y5,w,'omega');
         direction_der(x3,y3,x6,y6,w,'x1') direction_der(x3,y3,x6,y6,w,'y1') direction_der(x3,y3,x6,y6,w,'omega')];    

     %Normal matrix
     N = A'*P*A; 
     
     %Vector of right hand side of normal equations
     n = A'*P*l; 
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx = inv(N);

     %Solution of the normal equations
     x_hat = Q_xx*n;
       
     %Update
     X_hat=X_0+x_hat;
     X_0 = X_hat;

     x3 = X_0(1);
     y3 = X_0(2);
     w = X_0(3);

    
     %Check 1
     max_x_hat=max(abs(x_hat));
     
     %Vector of residuals
     v = A*x_hat-l;
 
     %Vector of adjusted observations
     L_hat = L+v;
    
     %Objective function
     vTPv = v'*P*v;
    
     %Functional relationships without the observations
     

     %Check 2
     L_0 = [direction_rad(x3,y3,x1,y1,w) direction_rad(x3,y3,x2,y2,w) direction_rad(x3,y3,x4,y4,w) direction_rad(x3,y3,x5,y5,w) direction_rad(x3,y3,x6,y6,w)]';
     Check2 = max(abs(L_hat-L_0));    
     
     %Update number of iterations
     iteration = iteration+1;
  
end

if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

X_hat

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2 * Q_xx; 

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));

%Cofactor matrix of adjusted observations
Q_LL_hat = A * Q_xx * A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2 * Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));

%Cofactor matrix of the residuals
Q_vv = Q_LL - Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));
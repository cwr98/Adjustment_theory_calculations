
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

r31 = 206.9094;
r32 = 46.5027;
r34 = 84.6449;
r35 = 115.5251;
r36 = 155.5891;

a342 = r34-r32
a354 = r35-r34
a365 = r36-r35
a316 = r31-r36
%Vector of observations (distances)
L = [a342 a354 a365 a316]'.*pi/200; %rad

%Number of observations
no_n = length(L);

%Initial values for the unknowns  
x3 = 200;
y3 = 460;

%Vector of initial values for the unknowns
X_0 = [x3 y3]';

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
S_rr = diag(s_L.^2);

F = [0 -1 1 0 0;
    0 0 -1 1 0;
    0 0 0 -1 1;
    1 0 0 0 -1];
S_LL = F*S_rr*F';

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
     L_0 = [angle_rad(x3,y3,x4,y4,x2,y2) angle_rad(x3,y3,x5,y5,x4,y4) angle_rad(x3,y3,x6,y6,x5,y5) angle_rad(x3,y3,x1,y1,x6,y6)]';
     
     %Vector of reduced observations
     l = L-L_0;
    
     %Design matrix with the elements from the Jacobian matrix J
     A = [angle_der(x3,y3,x4,y4,x2,y2,'xi') angle_der(x3,y3,x4,y4,x2,y2,'yi');
         angle_der(x3,y3,x5,y5,x4,y4,'xi') angle_der(x3,y3,x5,y5,x4,y4,'yi');
         angle_der(x3,y3,x6,y6,x5,y5,'xi') angle_der(x3,y3,x6,y6,x5,y5,'yi');
         angle_der(x3,y3,x1,y1,x6,y6,'xi') angle_der(x3,y3,x1,y1,x6,y6,'yi')];

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
     L_0 = [angle_rad(x3,y3,x4,y4,x2,y2) angle_rad(x3,y3,x5,y5,x4,y4) angle_rad(x3,y3,x6,y6,x5,y5) angle_rad(x3,y3,x1,y1,x6,y6)]';
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
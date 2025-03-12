
clc;
clear all;
close all;

%fixed values
xa = 0.0;
ya = 0.0;

xb = 0.0;
yb = 800.0;

%obs
a_cab = 29+(44/60)+(45/3600);
a_abc = 75+(57/60)+(53/3600);
a_bca = 74+(17/60)+(26/3600);

s_ac = 806.21;
s_bc = 412.32;

%Vector of observations
L = [a_cab*(pi/180) a_abc*(pi/180) a_bca*(pi/180) s_ac s_bc]'; %rad | m

%Number of observations
no_n = length(L);

%Initial values for the unknowns  
xc = 400;
yc = 700;

%Vector of initial values for the unknowns
X_0 = [xc yc]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n - no_u; 

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
ss = 10/100; %m
sr = (3/3600)*pi/180; %rad
s_L = [sr sr sr 1.5/100 1/100];
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
     L_0 = [angle_rad(xa,ya,xb,yb,xc,yc) angle_rad(xb,yb,xc,yc,xa,ya) angle_rad(xc,yc,xa,ya,xb,yb) distance_(xa,ya,xc,yc) distance_(xb,yb,xc,yc)]';

     %Vector of reduced observations
     l = L-L_0;
    
     %Design matrix with the elements from the Jacobian matrix J
     A = [angle_der(xa,ya,xb,yb,xc,yc,'xl') angle_der(xa,ya,xb,yb,xc,yc,'yl');
         angle_der(xb,yb,xc,yc,xa,ya,'xk') angle_der(xb,yb,xc,yc,xa,ya,'yk');
         angle_der(xc,yc,xa,ya,xb,yb,'xi') angle_der(xc,yc,xa,ya,xb,yb,'yi');
         distance_der(xa,ya,xc,yc,'x2') distance_der(xa,ya,xc,yc,'y2');
         distance_der(xb,yb,xc,yc,'x2') distance_der(xb,yb,xc,yc,'y2')];

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

     xc = X_0(1);
     yc = X_0(2);

    
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
     L_0 = [angle_rad(xa,ya,xb,yb,xc,yc) angle_rad(xb,yb,xc,yc,xa,ya) angle_rad(xc,yc,xa,ya,xb,yb) distance_(xa,ya,xc,yc) distance_(xb,yb,xc,yc)]';
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
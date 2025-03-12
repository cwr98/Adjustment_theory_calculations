
clc;
clear all;
close all;

%fixed values
x6 = 4968940.373;
y6 = 5317651.428;

x9 = 4970922.160;
y9 = 5324162.853;

%obs
r16 = 148.0875*pi/200;
r15 = 228.9044*pi/200;
r61 = 248.0883*pi/200;
r69 = 81.1917*pi/200;
r915 = 207.9027*pi/200;
r91 = 248.4428*pi/200;
r96 = 261.1921*pi/200;
r151 = 358.9060*pi/200;
r159 = 57.9014*pi/200;

s61 = 4307.851;
s91 = 10759.852;
s96 = 6806.332;
s151 = 6399.069;
s159 = 8751.757;

%Vector of observations
L = [r16 r15 r61 r69 r915 r91 r96 r151 r159 s61 s91 s96 s151 s159]'; %rad | m

%Number of observations
no_n = length(L);

%Initial values for the unknowns  
x1 = 4965804.18;
y1 = 5314698.13;
x15 = 4962997.53;
y15 = 5320448.85;
w1=1;
w6=1;
w9=1;
w15=1;

%Vector of initial values for the unknowns
X_0 = [x1 y1 x15 y15 w1 w6 w9 w15]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n - no_u; 

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
ss = 10/100; %m
sr = (1/1000)*pi/200; %rad
s_L = [sr sr sr sr sr sr sr sr sr ss ss ss ss ss];
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
     L_0 = [direction_rad(x1,y1,x6,y6,w1) direction_rad(x1,y1,x15,y15,w1) direction_rad(x6,y6,x1,y1,w6) direction_rad(x6,y6,x9,y9,w6) direction_rad(x9,y9,x15,y15,w9) direction_rad(x9,y9,x1,y1,w9) direction_rad(x9,y9,x6,y6,w9) direction_rad(x15,y15,x1,y1,w15) direction_rad(x15,y15,x9,y9,w15) distance_(x6,y6,x1,y1) distance_(x9,y9,x1,y1) distance_(x9,y9,x6,y6) distance_(x15,y15,x1,y1) distance_(x15,y15,x9,y9)]';
     
     %Vector of reduced observations
     l = L-L_0;
    
     %Design matrix with the elements from the Jacobian matrix J
     A = [direction_der(x1,y1,x6,y6,w1,'x1') direction_der(x1,y1,x6,y6,w1,'y1') 0 0 -1 0 0 0;
         direction_der(x1,y1,x15,y15,w1,'x1') direction_der(x1,y1,x15,y15,w1,'y1') direction_der(x1,y1,x15,y15,w1,'x2') direction_der(x1,y1,x15,y15,w1,'y2') -1 0 0 0;
         direction_der(x6,y6,x1,y1,w6,'x2') direction_der(x6,y6,x1,y1,w6,'y2') 0 0 0 -1 0 0;
         0 0 0 0 0 -1 0 0;
         0 0 direction_der(x9,y9,x15,y15,w9,'x2') direction_der(x9,y9,x15,y15,w9,'y2') 0 0 -1 0;
         direction_der(x9,y9,x1,y1,w9,'x2') direction_der(x9,y9,x1,y1,w9,'y2') 0 0 0 0 -1 0;
         0 0 0 0 0 0 -1 0;
         direction_der(x15,y15,x1,y1,w15,'x2') direction_der(x15,y15,x1,y1,w15,'y2') direction_der(x15,y15,x1,y1,w15,'x1') direction_der(x15,y15,x1,y1,w15,'y1') 0 0 0 -1;
         0 0 direction_der(x15,y15,x9,y9,w15,'x1') direction_der(x15,y15,x9,y9,w15,'y1') 0 0 0 -1;
         distance_der(x6,y6,x1,y1,'x2') distance_der(x6,y6,x1,y1,'y2') 0 0 0 0 0 0;
         distance_der(x9,y9,x1,y1,'x2') distance_der(x9,y9,x1,y1,'y2') 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0;
         distance_der(x15,y15,x1,y1,'x2') distance_der(x15,y15,x1,y1,'y2') distance_der(x15,y15,x1,y1,'x1') distance_der(x15,y15,x1,y1,'y1') 0 0 0 0;
         0 0 distance_der(x15,y15,x9,y9,'x1') distance_der(x15,y15,x9,y9,'y1') 0 0 0 0];

        

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

     x1 = X_0(1);
     y1 = X_0(2);
     x15 = X_0(3);
     y15 = X_0(4);
     w1=X_0(5);
     w6=X_0(6);
     w9=X_0(7);
     w15=X_0(8);

    
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
     L_0 = [direction_rad(x1,y1,x6,y6,w1) direction_rad(x1,y1,x15,y15,w1) direction_rad(x6,y6,x1,y1,w6) direction_rad(x6,y6,x9,y9,w6) direction_rad(x9,y9,x15,y15,w9) direction_rad(x9,y9,x1,y1,w9) direction_rad(x9,y9,x6,y6,w9) direction_rad(x15,y15,x1,y1,w15) direction_rad(x15,y15,x9,y9,w15) distance_(x6,y6,x1,y1) distance_(x9,y9,x1,y1) distance_(x9,y9,x6,y6) distance_(x15,y15,x1,y1) distance_(x15,y15,x9,y9)]';
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
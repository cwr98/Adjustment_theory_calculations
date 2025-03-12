
clc;
clear all;
close all;

%fixed values
x1 = 35.161;
y1 = 18.736;

x2 = 101.010;
y2 = 45.974;

x3 = 88.260;
y3 = 93.030;

x4 = 36.450;
y4 = 102.702;

x5 = 16.705;
y5 = 42.775;

a = 3.7/100;

%Vector of observations (distances)
L = [32.813 56.196 61.052 29.240 49.585 26.700 36.511 61.010]'; %m
%Lmod = L - ones(8,1)*a;

%Number of observations
no_n = length(L);

%Initial values for the unknowns  
x_100 = 40;
y_100 = 50;
x_200 = 60;
y_200 = 80;


%Vector of initial values for the unknowns
X_0 = [x_100 y_100 x_200 y_200]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n - no_u; 

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_L = [5 5 5 5 5 5 5 5]./1000;
S_LL = diag(s_L.^2);

%Theoretical standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_LL = (1/sigma_0^2)*S_LL;
%Q_LL =eye(9);

%Weight matrix
P = inv(Q_LL)
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
     L_0 = [distance_(x_100,y_100,x1,y1)+a distance_(x_100,y_100,x2,y2)+a distance_(x_100,y_100,x3,y3)+a distance_(x_100,y_100,x5,y5)+a distance_(x_200,y_200,x2,y2)+a distance_(x_200,y_200,x3,y3)+a distance_(x_200,y_200,x4,y4)+a distance_(x_200,y_200,x5,y5)+a]';
     
     %Vector of reduced observations
     l = L-L_0;

     
    
     %Design matrix with the elements from the Jacobian matrix J
     A = [distance_der(x_100,y_100,x1,y1,'x1') distance_der(x_100,y_100,x1,y1,'y1') 0 0;
         distance_der(x_100,y_100,x2,y2,'x1') distance_der(x_100,y_100,x2,y2,'y1') 0 0;
         distance_der(x_100,y_100,x3,y3,'x1') distance_der(x_100,y_100,x3,y3,'y1') 0 0;
         distance_der(x_100,y_100,x5,y5,'x1') distance_der(x_100,y_100,x5,y5,'y1') 0 0;
         0 0 distance_der(x_200,y_200,x2,y2,'x1') distance_der(x_200,y_200,x2,y2,'y1');
         0 0 distance_der(x_200,y_200,x3,y3,'x1') distance_der(x_200,y_200,x3,y3,'y1');
         0 0 distance_der(x_200,y_200,x4,y4,'x1') distance_der(x_200,y_200,x4,y4,'y1');
         0 0 distance_der(x_200,y_200,x5,y5,'x1') distance_der(x_200,y_200,x5,y5,'y1');];
    

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

     x_100 = X_0(1);
     y_100 = X_0(2);
     x_200 = X_0(3);
     y_200 = X_0(4);

    
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
     L_0 = [distance_(x_100,y_100,x1,y1)+a distance_(x_100,y_100,x2,y2)+a distance_(x_100,y_100,x3,y3)+a distance_(x_100,y_100,x5,y5)+a distance_(x_200,y_200,x2,y2)+a distance_(x_200,y_200,x3,y3)+a distance_(x_200,y_200,x4,y4)+a distance_(x_200,y_200,x5,y5)+a]';
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
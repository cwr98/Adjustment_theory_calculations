%--------------------------------------------------------------------------
%   
%          ADJUSTMENT CALCULATION I
%   Exercise 13: Adjustment Calculation - part VIII  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 15, 2018
%   Last changes   : October 15, 2018
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
%Load all files
dist = load('Distances.txt');
dir = load('Directions.txt');
fixedpoint = load('FixedPoints.txt');
newpoint = load('NewPoints.txt');

%Vector of observations
L = [dist(:,3); dir(:,3)*pi/200];   %[gon]->[rad]

%Gauss-Krueger coordinates for control points [m]

%Convert numbers to character array (num2str), 
%Execute MATLAB expression in text(eval)
for i=1:size(fixedpoint,1)
    eval(['y' num2str(fixedpoint(i,1)) '=' num2str(fixedpoint(i,2)) ';']);
    eval(['x' num2str(fixedpoint(i,1)) '=' num2str(fixedpoint(i,3)) ';']);
end
% y6 = fixedpoint(1,2);
% x6 = fixedpoint(1,3);
% y9 = fixedpoint(2,2);
% x9 = fixedpoint(2,3);

%Gauss-Krueger coordinates for new points [m]

%Convert numbers to character array (num2str), 
%Execute MATLAB expression in text(eval)
for i=1:size(newpoint,1)
    eval(['y' num2str(newpoint(i,1)) '=' num2str(newpoint(i,2)) ';']);
    eval(['x' num2str(newpoint(i,1)) '=' num2str(newpoint(i,3)) ';']);
end
% y1 = newpoint(1,2);
% x1 = newpoint(1,3);
% y15 = newpoint(2,2);
% x15 = newpoint(2,3);

%initial values for orientation unknowns
w1 = 0;
w6 = 0;
w9 = 0;
w15 = 0;

%Initial values for unknowns
X_0 = [x1 y1 x15 y15 w1 w6 w9 w15]';

%Number of observations
no_n = length(L);

%Number of unknowns
no_u = length(X_0);
 
%Redundancy
r = no_n-no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_dist = 0.1;           %[m]
s_dir = 0.001*pi/200;   %[gon]->[rad]

s_LL = [s_dist^2*ones(length(dist),1); s_dir^2*ones(length(dir),1)];
S_LL = diag(s_LL);

%Theoretical standard deviation
sigma_0 = 1;

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off condition
epsilon = 10^-5;
delta = 10^-9;
max_x_hat = Inf;
Check2 = Inf;

%Number of iterations
iteration = 0;

%Initialising A
A = zeros(no_n,no_u);

%Iteration
while max_x_hat>epsilon || Check2>delta

    %Vector of reduced distances
    L_0(1) = distance(y6,x6,y1,x1);
    L_0(2) = distance(y9,x9,y1,x1);
    L_0(3) = distance(y9,x9,y6,x6);
    L_0(4) = distance(y15,x15,y1,x1);
    L_0(5) = distance(y15,x15,y9,x9);
    
    %Vector of reduced directions
    L_0(6) = direction(y1,x1,y6,x6,w1);
    L_0(7) = direction(y1,x1,y15,x15,w1);
    L_0(8) = direction(y6,x6,y1,x1,w6);
    L_0(9) = direction(y6,x6,y9,x9,w6);
    L_0(10) = direction(y9,x9,y15,x15,w9);
    L_0(11) = direction(y9,x9,y1,x1,w9);
    L_0(12) = direction(y9,x9,y6,x6,w9);
    L_0(13) = direction(y15,x15,y1,x1,w15);
    L_0(14) = direction(y15,x15,y9,x9,w15);
    
    %Vector of reduced observations
    l = L-L_0';

    %Design matrix with the elements from the Jacobian matrix J
    A(1,1) = ds_dx_to(y6,x6,y1,x1);
    A(1,2) = ds_dy_to(y6,x6,y1,x1);
    
    A(2,1) = ds_dx_to(y9,x9,y1,x1);
    A(2,2) = ds_dy_to(y9,x9,y1,x1);
    
    A(4,1) = ds_dx_to(y15,x15,y1,x1);
    A(4,2) = ds_dy_to(y15,x15,y1,x1);
    A(4,3) = ds_dx_from(y15,x15,y1,x1);
    A(4,4) = ds_dy_from(y15,x15,y1,x1);

    A(5,3) = ds_dx_from(y15,x15,y9,x9);
    A(5,4) = ds_dy_from(y15,x15,y9,x9);
    
    A(6,1) = dr_dx_from(y1,x1,y6,x6);
    A(6,2) = dr_dy_from(y1,x1,y6,x6);
    A(6,5) = -1;
    
    A(7,1) = dr_dx_from(y1,x1,y15,x15);
    A(7,2) = dr_dy_from(y1,x1,y15,x15);
    A(7,3) = dr_dx_to(y1,x1,y15,x15);
    A(7,4) = dr_dy_to(y1,x1,y15,x15);
    A(7,5) = -1;
    
    A(8,1) = dr_dx_to(y6,x6,y1,x1);
    A(8,2) = dr_dy_to(y6,x6,y1,x1);
    A(8,6) = -1;
    
    A(9,6) = -1;
    
    A(10,3) = dr_dx_to(y9,x9,y15,x15);
    A(10,4) = dr_dy_to(y9,x9,y15,x15);
    A(10,7) = -1;
    
    A(11,1) = dr_dx_to(y9,x9,y1,x1);
    A(11,2) = dr_dy_to(y9,x9,y1,x1);
    A(11,7) = -1;
    
    A(12,7) = -1;
    
    A(13,1) = dr_dx_to(y15,x15,y1,x1);
    A(13,2) = dr_dy_to(y15,x15,y1,x1);
    A(13,3) = dr_dx_from(y15,x15,y1,x1);
    A(13,4) = dr_dy_from(y15,x15,y1,x1);
    A(13,8) = -1;
    
    A(14,3) = dr_dx_from(y15,x15,y9,x9);
    A(14,4) = dr_dy_from(y15,x15,y9,x9);
    A(14,8) = -1;
    
    %Normal matrix
    N = A'*P*A;

    %Vector of right hand side of normal equations
    n = A'*P*l;

    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx = inv(N);
    
    %Solution of normal equations
    x_hat = Q_xx*n;
    
    %Adjusted unknowns
    X_hat = X_0+x_hat;
    
    %Update
    X_0 = X_hat;
    
    x1 = X_0(1);
    y1 = X_0(2);
    x15 = X_0(3);
    y15 = X_0(4);
    w1 = X_0(5);
    w6 = X_0(6);
    w9 = X_0(7);
    w15 = X_0(8);
    
    %Check 1
    max_x_hat = max(abs(x_hat));
    
     %Vector of residuals
     v = A*x_hat-l;
 
     %Vector of adjusted observations
     L_hat = L+v;
    
     %Objective function
     vTPv = v'*P*v;
    
     %Functional relationships without the observations
     phi_X_hat = [distance(y6,x6,y1,x1);
                  distance(y9,x9,y1,x1);
                  distance(y9,x9,y6,x6);
                  distance(y15,x15,y1,x1);
                  distance(y15,x15,y9,x9);
                  direction(y1,x1,y6,x6,w1);
                  direction(y1,x1,y15,x15,w1);
                  direction(y6,x6,y1,x1,w6);
                  direction(y6,x6,y9,x9,w6);
                  direction(y9,x9,y15,x15,w9);
                  direction(y9,x9,y1,x1,w9);
                  direction(y9,x9,y6,x6,w9);
                  direction(y15,x15,y1,x1,w15);
                  direction(y15,x15,y9,x9,w15)];

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

%Convert to [gon] and check the quadrants
gon = X_0(5:8,1)*200/pi;
gon(1) = gon(1)+400;
gon(2) = gon(2)+400;
gon(3) = gon(3);
gon(4) = gon(4)+400;

%Results for the unknowns
Res = [x1;y1;x15;y15;gon(1);gon(2);gon(3);gon(4)]    %[m]/[gon]

%Vector of residuals -> Convert to [gon]
v_gon = v(6:14,1)*200/pi;

%Vector of adjusted observations -> Convert to [gon]
L_hat_gon = L_hat(6:14)*200/pi;

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknows
s_X = sqrt(diag(S_XX_hat));
s_X_gon = s_X(5:8,1)*200/pi;      %Convert to [gon]

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));
s_L_hat_gon = s_L_hat(6:14)*200/pi;  %Convert to [gon]

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;              

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));
s_v_gon = s_v(6:14,1)*200/pi;   %Convert to [gon]



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
%   Task 2
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Observations and initial values for unknowns
%--------------------------------------------------------------------------
% Load all files
dist = load('Distances.txt'); 
dir = load('Directions.txt'); 
fixedpoint = load('FixedPoints.txt');      % y x 6 9
np = load('NewPoints.txt');                % y x 1 15
coordinates = [fixedpoint; np];
coordinates2 = [fixedpoint; np];
dist_nums = dist(:,1:2);
dir_nums = dir(:,1:2);  
u_nums = np(:,1);
w_nums = coordinates(:,1);
% Vector of observations
L = [dist(:,3); dir(:,3)*pi/200];

% Gauss-Krueger coordinates for control points [m]
% for i=1:size(fixedpoint,1)
%     eval(['y' num2str(fixedpoint(i,1)) '=' num2str(fixedpoint(i,2)) ';']);
%     eval(['x' num2str(fixedpoint(i,1)) '=' num2str(fixedpoint(i,3)) ';']);
% end
fixed = [fixedpoint(:, 2) fixedpoint(:, 3)];

% Gauss-Krueger coordinates for new points [m]
% for i=1:size(newpoint,1)
%     eval(['y' num2str(newpoint(i,1)) '=' num2str(newpoint(i,2)) ';']);
%     eval(['x' num2str(newpoint(i,1)) '=' num2str(newpoint(i,3)) ';']);
% end

% Initial values for orientation unknowns
initial_w = zeros(length(coordinates),1); % w6 w9 w1 w15

% Initial values for unknowns
X_0 = [np(:,3); np(:,2); initial_w];    
% X_0  = x100 x101 x102 x103 y100 y101 y102 y103 w1000 w2000 w3000 w100 w101 w102 w103

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

% VC Matrix of the observations
s_dir = ones(length(dist),1)*0.001*pi/200;      
s_dist = ones(length(dir),1)*0.001;             
S_dist = diag(s_dist.^2);
S_dir = diag(s_dir.^2);
S_LL = blkdiag(S_dist, S_dir);

% Theoretical standard deviation
sigma_0 = 1e-4;

% Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

% Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
% break-off condition
epsilon = 1e-5;
delta = 1e-9;
max_x_hat = Inf;
Check2=Inf;

% Number of iterations
iteration = 0;

% Initialising A
A = zeros(no_n,no_u);

% Iteration
while max_x_hat>epsilon || Check2>delta

	% Distances
	L_0 = create_L_0(dist_nums, dir_nums, coordinates, initial_w);
	
	% Directions
	%L_0=

    
    % Vector of reduced observations
    l = L-L_0';

    % Design matrix with the elements from the Jacobian matrix J
	A = create_design_matrix(u_nums, w_nums, coordinates, dist_nums, dir_nums, no_n,no_u, initial_w);
  
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
    X_0 = X_hat;
    coordinates(length(fixedpoint)+1:length(coordinates),3) = X_hat(1:length(u_nums));
    coordinates(length(fixedpoint)+1:length(coordinates),2) = X_hat(length(u_nums) +1:2*length(u_nums));
    initial_w = X_hat(2*length(u_nums)+1:length(X_hat));
    
    % Check 1
    max_x_hat = max(abs(x_hat));
    
    % Vector of residuals
    v = A*x_hat-l;
 
    % Vector of adjusted observations
    L_hat = L+v;
    
    % Function
    vTPv = v'*P*v;
    
    % Functional relationships 
    phi_X_hat = create_L_0(dist_nums, dir_nums, coordinates, initial_w);
	 
    % Check 2
    Check2_full = L_hat - phi_X_hat;
    Check2 = max(abs(L_hat - phi_X_hat)); 
    
    % Update number of iterations
    iteration = iteration+1;

end

if Check2<=delta
    disp('Everything is fine!')
 else
    disp('Something is wrong.')
end

% Convert to [gon]
 X_hat(2*length(u_nums)+1:length(X_hat)) = X_hat(2*length(u_nums)+1:length(X_hat))*200/pi;
 for i = (2*length(u_nums)+1):length(X_hat)
 if X_hat(i) < 0
        X_hat(i) = X_hat(i) + 400;
 end
 end




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

% Results
table(X_hat, s_X, 'RowNames', {'x100' 'x101' 'x102' 'x103' 'y100' 'y101' 'y102' 
'y103' 'w1000' 'w2000' 'w3000' 'w100' 'w101' 'w102' 'w103'})

table(L, v, L_hat, s_v, s_L_hat)

function A = create_design_matrix(u_nums, w_nums, coordinates, dist_nums, dir_nums, no_n, no_u, initial_w)

A = zeros(no_n,no_u);

C_dist = x_y_dist(dist_nums, coordinates);
C_dir = x_y_w_dir(dir_nums, coordinates, initial_w);
for i = 1:length(dist_nums) % distance Jacobian matrix
    for j = 1:length(u_nums)
        if dist_nums(i, 1) == u_nums(j)
            A_dist(i,j) = ds_dx_from(C_dist(i,3), C_dist(i,1), C_dist(i,4), C_dist(i,2));
            A_dist(i,j+length(u_nums)) = ds_dy_from(C_dist(i,3), C_dist(i,1), C_dist(i,4), C_dist(i,2));
        elseif dist_nums(i, 2) == u_nums(j)
            A_dist(i,j) = ds_dx_to(C_dist(i,3), C_dist(i,1), C_dist(i,4), C_dist(i,2));
        A_dist(i,j+length(u_nums)) = ds_dy_to(C_dist(i,3), C_dist(i,1), C_dist(i,4), C_dist(i,2));
        end
    end
end
for i = 1:length(dir_nums) % directions Jacobian matrix
    for j = 1:length(u_nums)
        if dir_nums(i, 1) == u_nums(j)
            A_dir(i,j) = dr_dx_from(C_dir(i,3), C_dir(i,1), C_dir(i,4), C_dir(i,2));
            A_dir(i,j+length(u_nums)) = dr_dy_from(C_dir(i,3), C_dir(i,1), C_dir(i,4), C_dir(i,2));
            elseif dir_nums(i, 2) == u_nums(j)
                A_dir(i,j) = dr_dx_to(C_dir(i,3), C_dir(i,1), C_dir(i,4), C_dir(i,2));
                A_dir(i,j+length(u_nums)) = dr_dy_to(C_dir(i,3), C_dir(i,1), C_dir(i,4), C_dir(i,2));
            end
        end
end
for i = 1:length(dir_nums) % orientation parameter w Jacobian matrix
    for j = 1:length(w_nums)
        if dir_nums(i, 1) == w_nums(j)
            A_w(i,j) = -1;
        end
    end
end

A(1:size(A_dist,1), 1:size(A_dist,2)) = A(1:size(A_dist,1), 1:size(A_dist,2)) + A_dist;
A(size(A_dist,1)+1:size(A_dist,1)+size(A_dir,1), 1:size(A_dir,2)) = A(size(A_dist,1)+1:size(A_dist,1)+size(A_dir,1),  1:size(A_dir,2):size(A_dir,2)) + A_dir;
A(size(A_dist,1)+1:size(A_dist,1)+size(A_w,1), size(A_dir,2)+1:size(A_dir,2)+size(A_w,2)) = A(size(A_dist,1)+1:size(A_dist,1)+size(A_w,1), size(A_dir,2)+1:size(A_dir,2)+size(A_w,2)) + A_w;

function C = x_y_dist(dist_nums, coordinates)

for i = 1:size(dist_nums, 1)
    for j = 1:size(dist_nums, 2)
        for k = 1:size(coordinates, 1)
            if dist_nums(i,j) == coordinates(k, 1)
                Y(i,j) = coordinates(k,2);
                X(i,j) = coordinates(k,3);
            end
        end
    end
end
C = [X Y];    
end
end
function C = x_y_w_dir(G, coordinates, w)
for i = 1:size(G, 1)
    for j = 1:size(G, 2)
        for k = 1:size(coordinates, 1)
            if G(i,j) == coordinates(k,1) 
                Y(i,j) = coordinates(k,2);
                X(i,j) = coordinates(k,3); 
            end
            if G(i,1) == coordinates(k,1) 
                W(i,1) = w(k);
            end
        end
    end
end
C = [X Y W];    
end

function L_0 = create_L_0(dist_nums, dir_nums, coordinates, initial_w)

% Distances
L_0_dist = zeros(length(dist_nums),1);
C_dist = x_y_dist(dist_nums, coordinates);
for i = 1:length(dist_nums)
    L_0_dist(i) = distance(C_dist(i,3),C_dist(i,1),C_dist(i,4),C_dist(i,2));
end

% Directions
L_0_dir = zeros(length(dir_nums),1);

C_dir = x_y_w_dir(dir_nums, coordinates, initial_w);

for i = 1:length(dir_nums)
    L_0_dir(i) = direction(C_dir(i,3), C_dir(i,1), C_dir(i,4), C_dir(i,2), C_dir(i,5));
end

L_0 = [L_0_dist; L_0_dir];
end

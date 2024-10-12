%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%    Exercise 1: Introduction to Octave/Matlab   
% 
%   Author         : Anastasia Pasioti
%   Version        : October 02, 2020
%   Last changes   : October 02, 2021
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Task 1
%--------------------------------------------------------------------------
x = 4^2;
y = sqrt(4); 
z = 4^-2;
t = factorial(4);

%--------------------------------------------------------------------------
%   Task 2 - right-angled triangle
%--------------------------------------------------------------------------
%Using the Law of Sines
a = 4;                  %[m]
alpha = 53.1301;       %[degrees]
alpha_rad = (alpha/180) * pi;
c = a / cos(alpha_rad);                  %[m] 
disp('C:');
disp(c);

% Using the Pythagorean Theorem   
b = 3;                  %[m]
c_2 = sqrt(a^2+b^2);                %[m]
disp(c_2);

%--------------------------------------------------------------------------
%   Task 3 - General triangle
%--------------------------------------------------------------------------
%Using the law of cosines
a = 3;
b = 2;
gamma = 55;
c = sqrt(a^2+b^2-2*a*b*cosd(gamma));
disp(c)

%Using the law of sines
a = 4;
b = 3;
alpha = 30;
sin_beta = (b*sind(alpha))/a;
disp(sin_beta);
beta= asind(sin_beta);
disp(beta);

%--------------------------------------------------------------------------
%   Task 4 - Using matrices
%--------------------------------------------------------------------------
A = [1,2,3;4,5,6];
B = [4,6,5;8,2,3];
B_t = B'; 

%Calculate the product
C = A * B_t; 
disp('C product:');
disp(C);
%Load file and assign matrix to N
N = load('Matrix.txt');
%Calculate determinant, rank, inverse, pseudo-inverse, eigenvector and
%eigenvalues of matrix N
det_N = det(N);
%if you can create row from other row then it reduces ranks
rank_N = rank(N);
%if you multiply the inv with the orginal you get the identity matrix 
inv_N = inv(N);
%if we dint have inv we use this to get as close as possible 
pinv_N = pinv(N);
%eigenvectors = matrix * vector, gives vector that is scaled 
%eigenvalues = value of the scale 
[eigenvectors, eigenvalues] = eig(N);

disp('Determinant:');
disp(det_N);

disp('Rank:');
disp(rank_N);

disp('Inverse:');
disp(inv_N);

disp('Pseudo-inverse:');
disp(pinv_N);

disp('Eigenvalues:');
% Display eigenvalues as a vector
disp(diag(eigenvalues)); 
disp('Eigenvectors:');
disp(eigenvectors);







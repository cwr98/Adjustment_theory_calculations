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
clear variables;
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
a = 4;                   %[m]
alpha = 53.1301*pi/180;               %[degrees]
c=a/sin(alpha);                  %[m] 

% Using the Pythagorean Theorem   
b =3;                   %[m]
c_2 =sqrt(a^2+b^2);                 %[m]

%--------------------------------------------------------------------------
%   Task 3 - General triangle
%--------------------------------------------------------------------------
%Using the law of cosines
a=3;
b=3;
gamma=55*pi/180;
c=(a^2+b^2-2*a*b*cos(gamma));

%Using the law of sines
a=4;
b=3;
alpha=30*pi/180;
beta=asin(b/a*sin(alpha));
%always convert to gome

%--------------------------------------------------------------------------
%   Task 4 - Using matrices
%--------------------------------------------------------------------------
A = [1 2 3; 4 5 6];
B = [4 6 5; 8 2 3];

%Calculate the product
C=A*B';

%Load file and assign matrix to N
N = load("Matrix.txt");

%Calculate determinant, rank, inverse, pseudo-inverse, eigenvector and
%eigenvalues of matrix N
%3 blue one brow utube
det_N = det(N);
rank_N = rank(N);

%inv_N=inv(N);
inv_N=N^-1;
pseudoinv_N=pinv(N);

[V,D]=eig(N)








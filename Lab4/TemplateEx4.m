%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%   Exercise 4: Propagation of observation errors - part II  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 05, 2018
%   Last changes   : November 23, 2022
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Task 1
%--------------------------------------------------------------------------
disp('Task 1')

%Given
a = 3;           %[m]
b = 4;           %[m]
c = 5;           %[m]

s_abc = 0.02;    %[m]

%Semiperimeter of the triangle
s = (a+b+c)/2    %[m]
 
%Area of the triangle
A = sqrt(s*(s-a)*(s-b)*(s-c))        %[m^2]

%Functional relationship 


%Design matrices




%F = 

%Stochastic model
%S_LL = 

%VC propagation
%S_A = 

%Standard deviation
%s_a = 


%--------------------------------------------------------------------------
%   Task 2
%--------------------------------------------------------------------------
disp('Task 2')

%Given
s1 = 8;           %[m]
s2 = 6;           %[m]
t1 = 0*pi/200;    %[gon]->[rad]
t2 = 100*pi/200;  %[gon]->[rad]

s_s = 0.001;      %[m]
s_t = 0.1*pi/200; %[gon]->[rad]

%Functional relationships


%Design matrices




%F = 

%Stochastic model
%S_LL = 

%VC propagation
%S_XX = 

%Standard deviation
%s_x =      %[m]

%--------------------------------------------------------------------------
%   Task 3
%--------------------------------------------------------------------------
disp('Task 3')

%Given
r1 = 10;     %[m]
r2 = 14;     %[m]
H1 = 8;      %[m]
H2 = 12;     %[m]

srH = 0.01;  %[m]

%Functional relationships


%Amount of water at each water tank
%V1s = 
%V2s = 

%Design matrices


%F = 

%Stochastic model
%S_LL = 

%VC propagation
%S_xx = 

%Standard deviation
%s_x =






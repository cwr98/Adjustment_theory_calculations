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
s = (a+b+c)/2;    %[m]
 
%Area of the triangle
A = sqrt(s*(s-a)*(s-b)*(s-c));        %[m^2]

%Functional relationship 
ap = s - a;
bp = s - b;
cp = s - c;
A_2 = s*ap*bp*cp;
%Design matrices
F1 = [1 0 0; 0 1 0; 0 0 1; 1/2 1/2 1/2];
F2 = [-1 0 0 1; 0 -1 0 1; 0 0 -1 1; 0 0 0 1];
F3 = [s*bp*cp s*ap*cp s*ap*bp ap*bp*cp];
F4 = 1/2*1/sqrt(A_2);
F = F4*F3*F2*F1;

%Stochastic model
S_LL = s_abc^2*eye(3);

%VC propagation
S_A =  F*S_LL*F';

%Standard deviation
 s_a = sqrt(S_A);


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
x1 = s1 * cos(t1);
y1 = s1 * sin(t1);
x2 = s2 * cos(t2);
y2 = s2 * sin(t2);

dX = x2 - x1;
dY = y2 - y1;

d_2= dX^2 + dY^2;

d = sqrt(d_2);

disp(['Distance between two points: ' num2str(d) ' m'])

%Design matrices
 F1=[cos(t1) -s1*sin(t1) 0 0;
    sin(t1) s1*cos(t1) 0 0;
    0 0 cos(t2) -s2*sin(t2);
    0 0 sin(t2) s2*cos(t2)];
F2=[-1 0 1 0;
    0 -1 0 1];
F3=[2*dX 2*dY];

F4=[1/(2*sqrt(d_2))];

F = F4*F3*F2*F1;

%Stochastic model
S_LL=[s_s^2 0 0 0;
    0 s_t^2 0 0;
    0 0 s_s^2 0;
    0 0 0 s_t^2];

%VC propagation
S_XX = F*S_LL*F';

%Standard deviation
s_x = sqrt(S_XX);      
%[m]
disp(['Standard deviation of distance: ' num2str(s_x) ' m'])

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
a = r1^2*H1+r2^2*H2;
b = r1^2+r2^2;
H = a / b;

%Amount of water at each water tank
V1s = pi*r1^2*H;
V2s = pi*r2^2*H; 

%Design matrices
%columns: r1 r2 H1 H2 rows: r1 r2 a b
F1=[1 0 0 0;   
   0 1 0 0;                       
   2*r1*H1 2*r2*H2 r1^2 r2^2;
   2*r1 2*r2 0 0];
%columns: r1 r2 a b rows: r1 r2 H
F2=[1 0 0 0;
   0 1 0 0;
   0 0 1/b -a/(b^2)];  
%columns: r1 r2 H rows: V1s V2s
F3=[2*pi*r1*H 0 pi*r1^2;
    0 2*pi*r2*H pi*r2^2];    

F = F3*F2*F1;

%Stochastic model
S_LL = srH^2*eye(4);

%VC propagation
S_XX = F*S_LL*F';

%Standard deviation
s_x = sqrt(S_XX);

disp(['Final volume of water in the 1st tank: ' num2str(V1s) ' m3';'Final volume of water in the 2nd tank: ' num2str(V2s) ' m3']);
disp(['SD of V1s: ' num2str(s_x(1,1)) ' m3']) 
disp(['SD of V2s: ' num2str(s_x(2,2)) ' m3'])







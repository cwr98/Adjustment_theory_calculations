%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%   Exercise 3: Propagation of observation errors - part I  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 05, 2018
%   Last changes   : November 16, 2022
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Task 1
%--------------------------------------------------------------------------
disp('Task 1')
syms a b

%Area of the rectangle
A = a*b;
%differentiate A in respect to a
dA_da=diff(A,a); %=b
%differentiate A in respect to b
dA_db=diff(A,b); %=a
%Given
a_val =15;        
b_val =25;        
sa =.03;       
sb =.04;       

%replace symbols with numbers  
A = subs(A,[a,b],[a_val,b_val]);
fprintf('A (value of A): %f\n', double(A));
dA_da = subs(dA_da, [a,b], [a_val,b_val]);
fprintf('dA/da (partial derivative w.r.t a): %f\n', double(dA_da));
dA_db = subs(dA_db, [a,b], [a_val,b_val]);
fprintf('dA/db (partial derivative w.r.t b): %f\n', double(dA_db));

%formula 
SA=sqrt(dA_da^2 * sa^2 + dA_db^2 * sb^2);
SA=eval(SA);
fprintf('STD: %f\n', double(SA));

%--------------------------------------------------------------------------
%   Task 2
%--------------------------------------------------------------------------
disp('Task 2')
syms b c alpha

%Area of the triangle
A =(sin(alpha)*b*c)/2;

%Given
b_val=15;
c_val=25;
alpha_val=55/200*pi;
sb=0.03;
sc=0.04;
s_alpha=0.1/200*pi;

%differentiate
dA_db=diff(A,b); %=
dA_dc=diff(A,c); %=
dA_d_alpha=diff(A,alpha); %=

A=eval(subs(A,[b,c,alpha],[b_val,c_val,alpha_val]));
fprintf('Area of the triangle: %f\n', double(A));
dA_db=eval(subs(dA_db,[b,c,alpha],[b_val,c_val,alpha_val]));
fprintf('dA/db (partial derivative w.r.t b): %f\n', double(dA_db));
dA_dc=eval(subs(dA_dc,[b,c,alpha],[b_val,c_val,alpha_val]));
fprintf('dA/dc (partial derivative w.r.t c): %f\n', double(dA_dc));
dA_d_alpha=eval(subs(dA_d_alpha,[b,c,alpha],[b_val,c_val,alpha_val]));
fprintf('dA/dalpha (partial derivative w.r.t alpha): %f\n', double(dA_d_alpha));

%Standard deviation
SA=sqrt(dA_d_alpha^2*s_alpha^2+dA_db^2*sb^2+dA_dc^2*sc^2);
fprintf('STD: %f\n', double(SA));
%--------------------------------------------------------------------------
%   Task 3
%--------------------------------------------------------------------------
disp('Task 3')
syms r sr
%Area of a circle 
A = pi*r^2;
%Given
r_val=100;
sr_val=0.01;

%differentiate
dA_dr=diff(A,r); %=

A=eval(subs(A,r,r_val));
fprintf('Area of the circle: %f\n', double(A));
dA_dr=eval(subs(dA_dr,r,r_val));
fprintf('dA/dr (partial derivative w.r.t r): %f\n', double(dA_dr));
%Standard deviation
SA=sqrt(dA_dr^2*sr_val^2);
fprintf('STD: %f\n', double(SA));
%--------------------------------------------------------------------------
%   Task 4
%--------------------------------------------------------------------------
disp('Task 4')
syms c 
%Radius of the circle
r =c/(2*pi);
%Given
c_val=0.3;
sigma_c = 0.001;

%differentiate
dR_dc=diff(r,c); %=

r=eval(subs(r,c,c_val));
fprintf('Radius of the circle: %f\n', double(r));
dR_dc=eval(subs(dR_dc,c,c_val));
fprintf('dR/dc (partial derivative w.r.t c): %f\n', double(dR_dc));
%Standard deviation
SA=sqrt(dR_dc^2*sigma_c^2);
fprintf('STD: %f\n', double(SA));
%--------------------------------------------------------------------------
%   Task 5
%--------------------------------------------------------------------------
disp('Task 5')

%Given
dist_org = 20;
sigmia_org = 0.004;
dist_new= 100;

%Standard deviation
s100 = sqrt((dist_new/dist_org) * sigmia_org^2);
fprintf('STD: %f\n', double(s100));
%--------------------------------------------------------------------------
%   Task 6
%--------------------------------------------------------------------------
disp('Task 6')
syms t
grav=9.8;
%Height of the main building
H =grav*(t^2)/2;
%Given
t_0=0;
t_end=2.98;
st=0.1;

%differentiate
dH_dt=diff(H,t); 
H=eval(subs([H,t], t_end-t_0));
fprintf('Hight of building and time in seconds: %f\n', double(H));

dH_dt = subs(dH_dt, t, t_end-t_0);
fprintf('dH/dt (partial derivative w.r.t t_0 to t_end): %f\n', double(dH_dt));

%Standard deviation
SA=eval(sqrt(dH_dt^2*st^2));
fprintf('STD: %f\n', double(SA));
%The height of the building is: 43.514m
%The standard deviation of the height of the building: 0.29204 m

%--------------------------------------------------------------------------
%   Task 7
%--------------------------------------------------------------------------
disp('Task 7')

% velocity (v) of car at a given time (t),\
% depending on start velocity (start_v) and acceleration (a)
syms t start_v a
v = start_v + a*t;

% Given numeric values
% (avoid overwriting symbols by adding _val suffix)
start_v_val=15;
a_val=2;
end_pos=1000;
sigma_start_v=0.2;
sigma_a=0.1;

% position of the car (pos) at some time (t1)
syms t1
pos=int(v,t, 0,t1);

% find the t1 that results in the position of the car (pos) being its
% end position (end_pos)
% (store the solution in t1_sol, to avoid overwriting symbols)
t1_sol=solve(pos==end_pos, t1);

% we have two solutions, we pick one of them, the right one after printing it:
t1_sol=-(start_v - (start_v^2 + 2000*a)^(1/2))/a;

% differentiate the solution of t1 with respect to start_v and a
dt1_dstart_v = diff(t1_sol, start_v);
dt1_da = diff(t1_sol, a);

% standard deviation for the solution of t1:
sigma_t1 = sqrt(dt1_dstart_v^2 * sigma_start_v^2 + dt1_da^2 * sigma_a^2);

% substitute the symbols with the actual numbers:
t1_sol_val=subs(t1_sol, [start_v, a], [start_v_val, a_val]);
fprintf('How much time does the car need to travel 1 km from the starting point 0: %f\n', double(t1_sol_val));
sigma_t1_val = eval(subs(sigma_t1, [start_v, a], [start_v_val, a_val]));
fprintf('STD for the solution of t1: %f\n', double(sigma_t1_val));







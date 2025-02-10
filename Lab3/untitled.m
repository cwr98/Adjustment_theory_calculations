clc;
clear all;
close all;

% syms time 
% 
% tanks=(time^3)/17;
% 
% time_val = 50;
% 
% d_tanks_d_time=diff(tanks,time) %=
% 
% factories_at_50=eval(subs(d_tanks_d_time,[time],[time_val]))
% 
% factories_at_100=eval(subs(d_tanks_d_time,[time],[100]))

% syms x y 
% f=7*x^5+1/7*y^4
% d_x= diff(f, x)
% d_y= diff(f, y)

% syms x
% y=x.^3+x.^2
% assume(x, 'real')
% x_val=vpa(solve(y==3, x))
% y=subs(y, {x}, {x_val})

% syms t 
% car_speed=2*t
% a=int(car_speed, t, 0, 20)

% syms x y z 
% 
% f=(4*x-(1/2*y.^2))*z;
% 
% df_dx=diff(f,x)
% df_dy=diff(f,y)
% df_dz=diff(f,z)

% syms t 
% 
% a=-9.81;
% v_0=5;
% h_0=1;
% 
% V_time=v_0+(int(a,t, 0, t))
% 
% height=h_0+int(V_time, t, 0, t)
% assume(t, 'positive')
% time_end=eval(solve(height == 0, t))

syms x a b 
f=x.^2;
diff(f,x)
int(f,x)
int(f,x, a,b)
int(f,x, 0,x)



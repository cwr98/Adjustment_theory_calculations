clc;
clear all;
close all;

format longG


s1 = 81.617; %m
s2 = 143.524; %m
a=1.473; %m
b=1.698;%m

z1 = 85.318*pi/200; %rad
z2 = 67.924*pi/200;

ha = 35.945;

hc = s1*cos(z1)-b+ha
hb=s2*cos(z2)-b+hc

s_l = [3/1000 5/1000 2/1000 2/1000 (8/1000)*pi/200 (12/1000)*pi/200];
s_ll = diag(s_l.^2)

f1 = [0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 0 1;
    cos(z1) 0 1 -1 -s1*sin(z1) 0];

f2 = [cos(z2) 1 -1 -s2*sin(z2) 1];

f = f2*f1;

sxx = f*s_ll*f'

sqrt(sxx)

x2 = s2-s3*cos(a2);
y2 = s3*sin(a2);
x4 = s1*cos(a1);
y4 = -s1*sin(a1);

d = distance_(x2,y2,x4,y4)

f1= [0 1 -cos(a2) 0 s3*sin(a2);
    0 0 sin(a2) 0 s3*cos(a2);
    cos(a1) 0 0 -s1*sin(a1) 0;
    -sin(a1) 0 0 -s1*cos(a1) 0];

f2 = [distance_der(x2,y2,x4,y4,'x1') distance_der(x2,y2,x4,y4,'y1') distance_der(x2,y2,x4,y4,'x2') distance_der(x2,y2,x4,y4,'y2')];

s_l = [1.2/100 1.9/100 3.6/100 (1.5/1000)*pi/200 (4.1/1000)*pi/200];
s_ll = diag(s_l.^2)

f = f2*f1;

s_xx = f*s_ll*f';
sqrt(diag(s_xx))


t1 = 9.7 %sec
t2 = 23.1 %sec

x1 = s1*sin(a1);
y1 = s1*cos(a1);
x2 = s2*sin(a2);
y2 = s2*cos(a2);

s = [1.2/100 1.9/100 3.6/100 (1.5/1000)*pi/200 (4.1/1000)*pi/200 0.1 0.1];
s_ll = diag(s.^2);

f1 = [s1*cos(a1) 0 sin(a1) 0 0 0;
    -s1*sin(a1) 0 cos(a1) 0 0 0;
    0 s2*cos(a2) 0 sin(a2) 0 0;
    0 -s2*sin(a2) 0 cos(a2) 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];

d12 = distance(x1, y1, x2, y2)

f2 = [distance_der(x1, y1, x2, y2, 'x1') distance_der(x1, y1, x2, y2, 'y1') distance_der(x1, y1, x2, y2, 'x2') distance_der(x1, y1, x2, y2, 'y2') 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1]

syms x1_param x2_param y1_param y2_param;
f = sqrt((x2_param-x1_param)^2+(y2_param-y1_param)^2);
df_dx1 = diff(f,x1_param);
df_dx2 = diff(f,x2_param);
df_dy1 = diff(f,y1_param);
df_dy2 = diff(f,y2_param);

% Specify specific values for v and a

% Substitute the values into the derivatives
df_dx1_value = double(subs(df_dx1, [x1_param, x2_param, y1_param, y2_param], [x1, x2, y1, y2]))
df_dy1_value = double(subs(df_dy1, [x1_param, x2_param, y1_param, y2_param], [x1, x2, y1, y2]))
df_dx2_value = double(subs(df_dx2, [x1_param, x2_param, y1_param, y2_param], [x1, x2, y1, y2]))
df_dy2_value = double(subs(df_dy2, [x1_param, x2_param, y1_param, y2_param], [x1, x2, y1, y2]))



%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%     Exercise 2: Fundamentals of statistics  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 04, 2018
%   Last changes   : October 31, 2022
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Task 3
%--------------------------------------------------------------------------
%f(x)=C-((x-2)/4)*C PDF described using unkown C 
syms x C b %variables  
f=C-((x-2)/4)*C; %line 1 in matlab 
Area = int(f,x,2,6); %using int to find area and using PDF,  
% x(postion in pdf), our range  
Cval=solve(Area == 1, C); %C=0.5 solving the C with Area=1 
 
f=Cval-((x-2)/4)*Cval; %true PDF with C 
%f(x)=3/4 - x/8 
 
F = int(f,x,2,b) %finding the distribution function 
%F(b)=-((b - 2)*(b - 10))/16 
 
E = eval(int(f*x, x,2,6)) %E=expectation/arithmetic mean or n-intfinty 
 
STD = eval(sqrt(int((x-E)^2*f, x,2,6)))%finding standard devtiation slide 33 
 
%question 1, 0 (expailn) 
%question 2, 0 (expailn) 
 
prob4_5=eval(int(f, x,4,5)) % What is the probability for a  
% realisation 𝑥 of a random variable 







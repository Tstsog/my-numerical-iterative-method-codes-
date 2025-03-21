% This Matlab code solves two-variable non-linear equation f(x) = 0 using Newton-Raphson iterative method [1].
%
% Ref. [1] S. Chapra, "Applied numerical methods with MATLAB", Mc Craw Hill, Singapore (2008).
%
% The system of non-linear equations: f1 = 5*x1^2 + x1*x2^2 + sin(2*x2)^2 - 2 = 0,  
%                                     f2 = exp(2*x1 - x2) + 4*x2 - 3 = 0.
%
% The Newton-Raphson iterative scheme: x^(k+1) = x^(k) - [Jacobian(x^(k)]^(-1)*f(x^(k)), where Jacobian, J(x) = df/dx, 
% and 'k' defines the k-th iteration.
%
% The first order derivative is taken with finite difference scheme.  
%
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% March 21, 2025 & University of North Dakota 
%
function [] = two_variables_newton_1
clc; clear two_variables_newton_1
%
format long
%
x1 = 1.0; tol = 1e-6;
x2 = 1.0;
dx = 0.001; iter_max = 10.;
%

for iter = 1:iter_max
    %
    f1 = function_f1(x1,x2);
    f2 = function_f2(x1,x2);
    f12_val = [f1;
               f2];
    x_val = [x1;
             x2];
    %
    [iter, x1, x2, f1, f2]
    %
    Jacobian_11 = (function_f1(x1+dx,x2) - function_f1(x1,x2))/dx; 
    Jacobian_12 = (function_f1(x1,x2+dx) - function_f1(x1,x2))/dx; 
    Jacobian_21 = (function_f2(x1+dx,x2) - function_f2(x1,x2))/dx; 
    Jacobian_22 = (function_f2(x1,x2+dx) - function_f2(x1,x2))/dx; 
    %
    Jab_mat = [Jacobian_11, Jacobian_12;
               Jacobian_21, Jacobian_22 ];
    %
    x_val = x_val - Jab_mat\f12_val; % x_n+1 = x_n - J^(-1)_n * f_n
    x1 = x_val(1);
    x2 = x_val(2);
    %
    if ((abs(f2)) <= tol)
        break;
    end
%
end

%%%
% [     iter,               x1,                  x2,              f1,                    f2]
%[1.000000000000000   1.000000000000000   1.000000000000000   4.826821810431806   3.718281828459045
% 2.000000000000000   0.617610659174055  -0.276086894227563   0.229432742792465   0.428309018814646
% 3.000000000000000   0.568220970392787  -0.313553610671231   0.014580426864662   0.008881107877005
% 4.000000000000000   0.567306329314917  -0.309434664323555   0.000034733532951   0.000074294682328
% 5.000000000000000   0.567297360887681  -0.309442277626739   0.000000059624354   0.000000092416060    ];


%%%
return
end
%
function f1 = function_f1(x1,x2)
%
f1 = 5*x1^2 + x1*x2^2 + sin(2*x2)^2 - 2;
return
end
%
function f2 = function_f2(x1,x2)
%
f2 = exp(2*x1 - x2) + 4*x2 - 3;
return
end

% This Matlab code solves three-variable non-linear equation f(x) = 0 using Newton-Raphson iterative method [1].
%
% Ref. [1] S. Chapra, "Applied numerical methods with MATLAB", Mc Craw Hill, Singapore (2008).
%
% The system of non-linear equations: 
%  f1 = 3*x1 - cos(x2*x3) - 1.5      = 0,  
%  f2 = 4*x1^2 - 625*x2^2 + 2*x3 - 1 = 0,
%  f3 = 20*x3 + exp(-x1*x2) + 9      = 0.
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
function [] = three_variables_newton_1
clc; clear two_variables_newton_1
%
format long 
%
x1 = 3.0; tol = 1e-6;
x2 = 2.0;
x3 = 1.5;
dx = 0.001; iter_max = 50.;
%

for iter = 1:iter_max
    f1 = function_f1(x1,x2,x3);
    f2 = function_f2(x1,x2,x3);
    f3 = function_f3(x1,x2,x3);
    f12_val = [f1;
               f2;
               f3];
    x_val = [x1;
             x2;
             x3];
    %
    [iter, x1, x2, x3, f1, f2, f3]
    %
    Jacobian_11 = (function_f1(x1+dx,x2,x3) - function_f1(x1,x2,x3))/dx; 
    Jacobian_12 = (function_f1(x1,x2+dx,x3) - function_f1(x1,x2,x3))/dx; 
    Jacobian_13 = (function_f1(x1,x2,x3+dx) - function_f1(x1,x2,x3))/dx; 
    %
    Jacobian_21 = (function_f2(x1+dx,x2,x3) - function_f2(x1,x2,x3))/dx; 
    Jacobian_22 = (function_f2(x1,x2+dx,x3) - function_f2(x1,x2,x3))/dx; 
    Jacobian_23 = (function_f2(x1,x2,x3+dx) - function_f2(x1,x2,x3))/dx; 
    %
    Jacobian_31 = (function_f3(x1+dx,x2,x3) - function_f3(x1,x2,x3))/dx; 
    Jacobian_32 = (function_f3(x1,x2+dx,x3) - function_f3(x1,x2,x3))/dx; 
    Jacobian_33 = (function_f3(x1,x2,x3+dx) - function_f3(x1,x2,x3))/dx;     
    %
    Jab_mat = [Jacobian_11, Jacobian_12, Jacobian_13;
               Jacobian_21, Jacobian_22, Jacobian_23;
               Jacobian_31, Jacobian_32, Jacobian_33];
    %
    x_val = x_val - Jab_mat\f12_val; % x_n+1 = x_n - J^(-1)_n * f_n
    x1 = x_val(1);
    x2 = x_val(2);
    x3 = x_val(3);
    %
    if ((abs(f2)) <= tol)
        break;
    end
%
end
%
% [iter, x1, x2, x3, f1, f2, f3]
%[1.0e+03 *
%0.001000000000000   0.003000000000000   0.002000000000000   0.001500000000000   0.008489992496600  -2.462000000000000 0.039002478752177
%1.0e+02 *
%0.020000000000000   0.004232278774166   0.009891506672804  -0.004511373008092  -0.011323913795098  -6.126976888699748 0.006351975366136
%1.0e+02 *
%0.030000000000000   0.008291655963695   0.004949378908557  -0.004765738573647   0.000151864314117  -1.523052827478585 0.001319173339260
%4.000000000000000   0.832730490584170   0.249000756990531  -0.489872467065727   0.005621643509941 -37.956845267998141 0.015285376913214
%5.000000000000000   0.833230761651787   0.127275230439749  -0.494749084676833   0.001674196555021  -9.336769337884332 0.004398310727280
%6.000000000000000   0.833275603390139   0.068789462380265  -0.497159266510244   0.000411549873202  -2.174419442303382 0.001105968999600
%7.000000000000000   0.833281028365887   0.043660443229536  -0.498202780144133   0.000079645175667  -0.410372910724955 0.000216828591572
%8.000000000000000   0.833281582445461   0.036215067059786  -0.498512617497545   0.000007710371864  -0.039299378764358 0.000021091907181
%9.000000000000000   0.833281613514197   0.035357204301647  -0.498548337166759   0.000000196939038  -0.000996119667985 0.000000536886301
%10.000000000000000  0.833281613817609   0.035334938827777  -0.498549264373076   0.000000002828061  -0.000014225766961 0.000000007670721
%11.000000000000000  0.833281613816771   0.035334620651760  -0.498549277623129   0.000000000039556  -0.000000198923281 0.000000000107265    ];


%%%
return
end
%
function f1 = function_f1(x1,x2,x3)
%
f1 = 3*x1 - cos(x2*x3) - 1.5;
return
end
%
function f2 = function_f2(x1,x2,x3)
%
f2 = 4*x1^2 - 625*x2^2 + 2*x3 - 1;
return
end
%
function f3 = function_f3(x1,x2,x3)
%
f3 = 20*x3 + exp(-x1*x2) + 9;
return
end

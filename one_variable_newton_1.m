% This Matlab code solves one-variable non-linear equation f(x) = 0 using Newton-Raphson iterative method [1].
%
% Ref. [1] S. Chapra, "Applied numerical methods with MATLAB", Mc Craw Hill, Singapore (2008).
%
% The non-linear equation: f(x) = 2*sin(x) - x^2 = 0.  
% The Newton-Raphson iterative scheme: x^(k+1) = x^(k) - Jacobian^(-1)*f(x), where Jacobian, J(x) = df/dx, 
% and 'k' defines the k-th iteration.
%
% The first order derivative is taken with finite difference scheme.  
%
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% March 21, 2025 & University of North Dakota 
%
function [] = one_variable_newton_1
clc; clear one_variable_newton_1
%
format long
%
x = 2.0; tol = 1e-6;
dx = 0.001; iter_max = 10.;
%
for it = 1:iter_max
    %
    f = functon_one_varaible(x);
    %
    [it, x, f]
    % central difference scheme
    df = (functon_one_varaible(x+dx) - functon_one_varaible(x-dx))/(2*dx);    
    %
    % forward difference 
%    df = (functon_one_varaible(x+dx) - functon_one_varaible(x))/dx;    
    %
    Jacobian = f/df;
    x = x - Jacobian;
    %
    if (abs(f) <= tol)
        break; 
    end
    %
end
%
% [it,                      x,                f1 ]
%[1.000000000000000   2.000000000000000  -2.181405146348637
% 2.000000000000000   1.548577682454594  -0.398586486444075
% 3.000000000000000   1.418010103610696  -0.034050908835208
% 4.000000000000000   1.404559939988076  -0.000359580943717
% 5.000000000000000   1.404414840973710  -0.000000041825047];

%%%
return
end
%
function f = functon_one_varaible(x)
%
f = 2*sin(x) - x.^2;
%%%
return
end

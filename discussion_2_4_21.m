%% Problem 1
clear; clc

p = [3 -1 1];
q = [-2 0 1];
pq = q - p;

syms t;

r = p + pq*t;
disp(r);

%% Problem 2
% consider the parametric equations x = t and y = 4t^2 of a planar curve
clear; clc

syms t;

x = t;
y = 4*t^2;

% show that they represent a planar curve y = 4x^2
% just plug t = x into y = 4*t^2

% find the unit tangent vector
r = [x y];
rp = diff(r, t);
utv = rp/norm(rp);
disp('unit tangent vector:'); disp(utv);

% find the slope of the tangent line at any point
% slope is y/x for derivative of r
slope = rp(2)/rp(1);
disp('slope:'); disp(slope);

% find the parametric equations of the tangent to the curve at t = 1
p = double(subs(r, t, 1));
p_eqns = p + rp;

%% Problem 3
clear; clc

syms t;

r = [2*t, exp(2*t) 1];

% find integral
R = int(r);

% from 0 to 1
value = subs(R, t, 1) - subs(R, t, 0);
disp('value:'); disp(value);

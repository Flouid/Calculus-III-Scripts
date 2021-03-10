%% Problem 1
% conceptual

%% Problem 2
% Find the limit, if it exists. (If an answer does not exist, enter DNE.)
% lim (x, y) → (0, 0) (x^4 − 16y^2)/(x^2 + 8y^2)
clear; clc

syms x y;

% declare f(x, y)
f = (x^4 - 16*y^2)/(x^2 + 8*y^2);
% substitute y = 0 in to take the limit as we approach from the x direction
fx = subs(f, y, 0);
fxl = double(limit(fx, x, 0));
% substitute x = 0 in to take the limit as we approach from the y direction
fy = subs(f, x, 0);
fyl = double(limit(fy, y, 0));
% check if the limits are equal
if(fxl == fyl)
    disp('limit:'); disp(fxl);
else
    disp('DNE');
end

%% Problem 3
% Find the limit, if it exists. (If an answer does not exist, enter DNE.)
% lim (x, y) → (0, 0) x*y/(x^2 + y^2)^(1/2)
clear; clc

syms x y;

% declare f(x, y)
f = x*y/(x^2 + y^2)^(1/2);
% substitute y = 0 in to take the limit as we approach from the x direction
fx = subs(f, y, 0);
fxl = double(limit(fx, x, 0));
% substitute x = 0 in to take the limit as we approach from the y direction
fy = subs(f, x, 0);
fyl = double(limit(fy, y, 0));
% check if the limits are equal
if(fxl == fyl)
    disp('limit:'); disp(fxl);
else
    disp('DNE');
end

%% Problem 4
% Find h(x, y) = g(f(x, y)).
% g(t) = t^2 + t^(1/2),    f(x, y) = 7x + 5y − 35
% Find the set on which h is continuous.
clear; clc

syms x y t;

% declare g(t) and f(x, y)
g = t^2 + sqrt(t);
f = 7*x + 5*y - 35;
% substitute f(x, y) into g(t) to get h(x, y)
t = f;
h = subs(g);
% h is continuous where it is defined, since g and f are both continuous
% it is not defined where 7x + 5y - 35 < 0 because of the radical
ineq = f > 0;
% solve for y
domain = solve(ineq, y);

% output
disp('domain:'); disp(domain);

%% Problem 5
% Find h(x, y) = g(f(x, y)).
% g(t) = t + ln(t),    f(x, y) = (3 − x*y)/(2 + x^2*y^2)
clear; clc

syms x y t;

% declare g(t) and f(x, y)
g = t + log(t);
f = (3 - x*y)/(2 + x^2*y^2);
% substitute f(x, y) into g(t) to get h(x, y)
t = f;
h = subs(g);
% h is continuous where it is defined, since g and f are both continuous
% ln(x) <= 0 is not defined though. The denominator of f cannot be zero, so
% only the numerator matters. We wan't to find when the numerator is less
% than 3
[numerator, denominator] = numden(f);
ineq = numerator > 0;
% no further simplification needed

% output
disp('domain:'); disp(ineq);

%% Problem 6
% this one is just not meant for matlab

%% Problem 7
% Find the limit, if it exists. (If an answer does not exist, enter DNE.)
% lim (x, y) → (2, 1) 9x^3 − x^2*y^2
clear; clc

syms x y;

% declare f(x, y)
f = 9*x^3 - x^2*y^2;
% declare the given x and y
x = 2;
y = 1;
% for this limit, straight substitution will work
limit = double(subs(f));

% output
disp('limit:'); disp(limit);

%% Problem 8
% Find the limit, if it exists. (If an answer does not exist, enter DNE.)
% lim (x, y) → (−3, 3) e^(−x*y)*cos(x + y)
clear; clc

syms x y;

% declare f(x, y)
f = exp(-x*y)*cos(x + y);
% declare the given x and y
x = -3;
y = 3;
% for this limit, straight substitution will work
limit = double(subs(f)); % this is e^9

% output
disp('limit:'); disp('e^9');

%% Problem 9
% Find the limit, if it exists. (If an answer does not exist, enter DNE.)
% lim (x, y) → (0, 0) x^2*y*e^y/(x^4 + 8*y^2)
clear; clc

syms x y;

% declare f(x, y)
f = x^2*y*exp(y)/(x^4 + 8*y^2);
% substitute y = 0 in to take the limit as we approach from the x direction
fx = subs(f, y, 0);
fxl = double(limit(fx, x, 0));
% substitute x = 0 in to take the limit as we approach from the y direction
fy = subs(f, x, 0);
fyl = double(limit(fy, y, 0));
% substitute y = x^2 in to take the limit as we approach from a parabola
y = x^2;
fyx2 = subs(f);
fyx2l = subs(fyx2, x, 0);
% we have three limits, they all must equal
if((fxl == fyl) && (fxl == fyx2l))
    disp('limit:'); disp(fxl);
else
    disp('DNE');
end

%% Problem 10
% this problem is entirely visual, but trivial

%% Problem 11
% If f(x, y) = 4 − 2*x^2 − y^2, find fx(−4, 8) and fy(−4, 8) and interpret 
% these numbers as slopes.
clear; clc

syms x y;

% declare f(x, y)
f = 4 - 2*x^2 - y^2;
% take the partial derivatives with respect to x and y
fx = diff(f, x);
fy = diff(f, y);
% declare x and y as given
x = -4;
y = 8;
% substitute into the partial derivatives
fxv = double(subs(fx));
fyv = double(subs(fy));

% output
disp('fx(-4, 8):'); disp(fxv);
disp('fy(-4, 8):'); disp(fyv);

%% Problem 12
% Find fx(1, 0) and fy(1, 0) and interpret these numbers as slopes for the
% following equation. f(x, y) = (4 − x^2 − 2*y^2)^(1/2)
clear; clc

syms x y

% declare f(x, y)
f = (4 - x^2 - 2*y^2)^(1/2);
% take partial derivatives with respect to x and y
fx = diff(f, x);
fy = diff(f, y);
% declare x and y as given
x = 1;
y = 0;
% substitute into the partial derivatives
fxv = double(subs(fx));
fyv = double(subs(fy));

% output
disp('fx(1, 0):'); disp(subs(fxv));
disp('fy(1, 0):'); disp(fyv);

%% Problem 13
% Find the first partial derivatives of the function. 
% f(x, y) = x^8 + 9*x*y^7
clear; clc

syms x y;

% declare f(x, y)
f = x^8 + 9*x*y^7;
% take partial derivatives with respect to x and y
fx = diff(f, x);
fy = diff(f, y);

% output
disp('fx(x, y):'); disp(fx);
disp('fy(x, y):'); disp(fy);

%% Problem 14
% Find the first partial derivatives of the function. 
% f(x, t) = t^4*e^(−x)
clear; clc

syms x t;

% declare f(x, y)
f = t^4*exp(-x);
% take partial derivatives with respect to x and y
fx = diff(f, x);
ft = diff(f, t);

% output
disp('fx(x, y):'); disp(fx);
disp('fy(x, y):'); disp(ft);

%% Problem 15
% Find the first partial derivatives of the function.
% f(x, y) = x/y
clear; clc

syms x y;

% declare f(x, y)
f = x/y;
% take partial derivatives with respect to x and y
fx = diff(f, x);
fy = diff(f, y);

% output
disp('fx(x, y):'); disp(fx);
disp('fy(x, y):'); disp(fy);

%% Problem 16
% Find the first partial derivatives of the function.
% f(x, y) = (a*x + b*y)/(c*x + d*y)
clear; clc

syms x y a b c d;

% declare f(x, y)
f = (a*x + b*y)/(c*x + d*y);
% take partial derivatives with respect to x and y
fx = diff(f, x);
fy = diff(f, y);

% output
disp('fx(x, y):'); disp(fx);
disp('fy(x, y):'); disp(fy);

%% Problem 17
% Find the first partial derivatives of the function.
% w = ln(x + 9*y + 5*z)
clear; clc

syms x y z;

% declare w(x, y, z)
w = log(x + 9*y + 5*z);
% take partial derivatives with respect to x, y, and z
wx = diff(w, x);
wy = diff(w, y);
wz = diff(w, z);

% output
disp('wx(x, y, z);'); disp(wx);
disp('wy(x, y, z);'); disp(wy);
disp('wz(x, y, z);'); disp(wz);

%% Problem 18-19
% conceptual

%% Problem 20
% Find all the second partial derivatives.
% f(x, y) = x^6*y - 2*x^3*y^2
clear; clc

syms x y;

% declare f(x, y)
f = x^6*y - 2*x^3*y^2;
% take the first partial derivatives with respect to x and y
fx = diff(f, x);
fy = diff(f, y);
% take all of the second partial derivatives with respect to x and y
fxx = diff(fx, x);
fxy = diff(fx, y);
fyx = diff(fy, x);
fyy = diff(fy, y);

% output
disp('fxx(x, y):'); disp(fxx);
disp('fxy(x, y):'); disp(fxy);
disp('fyx(x, y):'); disp(fyx);
disp('fyy(x, y):'); disp(fyy);

%% Problem 21
% Find all the second partial derivatives.
% T = e^(−3*r)*cos(o)
clear; clc

syms r o;

% declare T(r, o)
T = exp(-3*r)*cos(o);
% take the first partial derivatives with respect to r and o
Tr = diff(T, r);
To = diff(T, o);
% take all of the second partial derivatives with respect to r and o
Trr = diff(Tr, r);
Tro = diff(Tr, o);
Tor = diff(To, r);
Too = diff(To, o);

% output
disp('Trr(r, o):'); disp(Trr);
disp('Tro(r, o):'); disp(Tro);
disp('Tor(r, o):'); disp(Tor);
disp('Too(r, o):'); disp(Too);
 
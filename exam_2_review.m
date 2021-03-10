% Louis Keith
% 2/27/21
% MA 231 Exam 2 Review
%% Problem 1
clear; clc
% Find the curvature of the following curve at t = 0.
% r(t) = <8*e^t, 2*t^2>

syms t;
r(t) = [8*exp(t) 2*t^2];
t_v = 0;

% MATLAB won't accept 2D vectors for cross products, so append a 0
r(t) = [r(t) 0];
% FORMULA: k = ||T'(t)||/||r'(t)|| = ||r'(t) x r''(t)||/||r'(t)||^3
r_p(t) = diff(r(t));            % find first derivative
r_dp(t) = diff(r_p(t));         % find second derivative
% apply second form of curvative formula
k(t) = norm(cross(r_p(t), r_dp(t)))/norm(r_p(t))^3;
% plug in value for t to get answers
k_0 = subs(k(t), t, t_v);
% output
disp('k(0):'); disp(k_0);

%% Problem 2
clear; clc
% Find the curvature of the curve y = cos(x) at x = π/6.

syms x t;
y = cos(x);
t_v = pi/6;

% parametrize the curve
r(t) = subs([t, y], x, t);
% COPY/PASTE the code from the previous problem

% MATLAB won't accept 2D vectors for cross products, so append a 0
r(t) = [r(t) 0];
% FORMULA: k = ||T'(t)||/||r'(t)|| = ||r'(t) x r''(t)||/||r'(t)||^3
r_p(t) = diff(r(t));            % find first derivative
r_dp(t) = diff(r_p(t));         % find second derivative
% apply second form of curvative formula
k(t) = norm(cross(r_p(t), r_dp(t)))/norm(r_p(t))^3;
% plug in value for t to get answers
k_0 = subs(k(t), t, t_v);
% output
disp('k(0):'); disp(k_0);

%% Problem 3
clear; clc
% Find the unit normal vector N(t) and the binormal vector vector B(t)
% to the curve vector r(t) = <5*cos(t), 5*sin(t), 5*ln(cos(t))>
% at (5, 0, 0).

syms t;
%r(t) = [5*cos(t) 5*sin(t) 5*log(cos(t))];
%r_v = [5, 0, 0];
r(t) = [t, 9*cos(t), 9*sin(t)];
r_v = [0, 9, 0];

% find t such that r(t) = r_v
t_v = solve(r(t) == r_v, t);
% FORMULA: T(t) = r'(t)/||r'(t)||
% FORMULA: N(t) = T'(t)/||T'(t)||
% FORMULA: B(t) = T(t) x N(t)
r_p(t) = diff(r(t));
T(t) = r_p(t)/norm(r_p(t));
T_p(t) = diff(T(t));
N(t) = T_p(t)/norm(T_p(t));
B(t) = cross(T(t), N(t));
% plug in value for t to get answers
N_t = subs(N(t), t, t_v);
B_t = subs(B(t), t, t_v);
% output
disp('N_t:'); disp(N_t);
disp('B_t:'); disp(B_t);

%% Problem 4
clear; clc
% Find the tangential and normal components of the acceleration of a 
% particle that moves along a path given by vector r(t) =  <t^2, e^t, 2*t>
% at t = 0.

syms t;
r(t) = [t^2 exp(t) 2*t];
t_v = 0;

% FORMULA: a_T = T(t) • r''(t) = |r''(t) • r'(t)|/|r'(t)|
% FORMULA: a_N = |T(t) x r''(t)| = |r'(t) x r''(t)||/|r'(t)|
r_p(t) = diff(r(t));        % the velocity vector
r_dp(t) = diff(r_p(t));     % the acceleration vector
a_T(t) = norm(dot(r_dp(t), r_p(t)))/norm(r_p(t));
a_N(t) = norm(cross(r_p(t), r_dp(t)))/norm(r_p(t));
% plug in value for t to get answers
a_T_v = subs(a_T(t), t, t_v);
a_N_v = subs(a_N(t), t, t_v);
% output
disp('a_T:'); disp(a_T_v);
disp('a_N:'); disp(a_N_v);

%% Problem 5
clear; clc
% This is a really easy domain recognition problem.

%% Problem 6
clear; clc
% Consider the function f(x,y) = (2*x^2 + 7*y^2)/(x^2 + y^2).
% Find the limit of f(x,y) when (x,y)→ (0,0) along each of the following 
% paths. 
% (a) x = 0
% (b) y = 0
% Does lim (x,y)→(0,0) f(x,y) exist?

syms x y;
fxy = (2*x^2 + 7*y^2)/(x^2 + y^2);
lim_target = [0 0];
xlim_target = lim_target(1);
ylim_target = lim_target(2);

% (a)
fy = subs(fxy, x, xlim_target);
lim1 = subs(fy, y, ylim_target);
% (b)
fx = subs(fxy, y, ylim_target);
lim2 = subs(fx, x, xlim_target);
% output
disp('(a):'); disp(lim1);
disp('(b):'); disp(lim2);
if (lim1 == lim2)
    fprintf('The limit is %d\n', lim1);
else
    disp('The limit does not exist');
end

%% Problem 7
clear; clc
% Consider the function f(x,y) = (9*x*y^3)/(x^4 + y^4).
% Find the limit of f(x,y) when (x,y)→ (0,0) along each of the following 
% paths. 
% (a) x = 0
% (b) y = 0
% Does lim (x,y)→(0,0) f(x,y) exist?

syms x y;
fxy = (9*x*y^3)/(x^4 + y^4);
lim_target = [0 x];
xlim_target = lim_target(1);
ylim_target = lim_target(2);

% (a)
fy = subs(fxy, x, xlim_target);
lim1 = subs(fy, y, ylim_target);
% (b)
fx = subs(fxy, y, ylim_target);
lim2 = subs(fx, x, xlim_target);
% output
disp('(a):'); disp(lim1);
disp('(b):'); disp(lim2);
if (lim1 == lim2)
    fprintf('The limit is %d\n', lim1);
else
    disp('The limit does not exist');
end

%% Problem 8
clear; clc
% Consider the function f(x,y) = (8*x^4*y^3)/(x^4 + y^4).

% SQUEEZE THEOREM PROBLEM
% x^4/(x^4 + y^4) is always between 0 and 1 according to arachchi
% meaning 8y^3 is always between 0 and 1.
% lim y->0 8y^3 = 0
% 0 <= lim <= 0

%% Problem 9
clear; clc
% Given z = x^2 + 7*x*y^3 + x - y + 2, find the partial derivatives
% ∂z/∂x and ∂z/∂y.

syms x y;
z = x^2 + 7*x*y^3 + x - y + 2;

% take partial derivatives
dzdx = diff(z, x);
dzdy = diff(z, y);
% output
disp('∂z/∂x:'); disp(dzdx);
disp('∂z/∂y:'); disp(dzdy);

%% Problem 10
clear; clc
% Given w = r^2*cos(5*t) + e^(r*sin(t)), find the partial derivatives
% ∂w/∂r and ∂w/∂t.

syms r t;
w = r^2*cos(5*t) + exp(r*sin(t));

% take partial derivatives
dwdr = diff(w, r);
dwdt = diff(w, t);
% output
disp('∂w/∂r:'); disp(dwdr);
disp('∂w/∂t:'); disp(dwdt);

%% Problem 11
clear; clc
% Given f(x, y) = (4*x*y)/(x^2 + y^4), find fx(1, 1) and fy(-1, 1).

syms x y;
f = (4*x*y)/(x^2 + y^4);
fx_v = [1 1];
fy_v = [-1 1];

% calculate fx(x, y) and fy(x, y)
fx = diff(f, x);
fy = diff(f, y);
% substitute values into fx(x, y) and fy(x, y) 
fx_subbed = subs(fx, [x y], fx_v);
fy_subbed = subs(fy, [x y], fy_v);
% output
disp('fx(1, 1):'); disp(fx_subbed);
disp('fy(-1, 1):'); disp(fy_subbed);

%% Problem 12
clear; clc
% Given f(x,y) = x^3*y + 9*e^(x*y), find fx, fxx, fy, fyy and fxy.

syms x y;
f = x^3*y + 9*exp(x*y);

% take derivatives
fx = diff(f, x);
fy = diff(f, y);
fxx = diff(fx, x);
fxy = diff(fx, y);
fyy = diff(fy, y);
% output
disp('fx(x):'); disp(fx);
disp('fxx(x):'); disp(fxx);
disp('fy(x):'); disp(fy);
disp('fyy(x):'); disp(fyy);
disp('fxy(x):'); disp(fxy);

%% Problem 13
clear; clc
% Find the equation of the tangent plane, in linear form, to the surface 
% f(x,y) = 8*x^2 − y^2 at the point (0, −2, −4).

syms x y z;
plane_eqn = z == 8*x^2 - y^2;
p = [0 -2 -4];

% EQUATION OF TANGENT PLANE:
% z − z0 = fx(x0, y0)*(x − x0) + fy(x0, y0)*(y − y0)
% to satisfy this, we need to take the partial derivatives with respect to
% x and y. Where z = f(x, y)
f = rhs(plane_eqn);
fx = diff(f, x);
fy = diff(f, y);
% get values by plugging in x and y from the specified point
fxv = double(subs(fx, [x y], [p(1) p(2)]));
fyv = double(subs(fy, [x y], [p(1) p(2)]));
% use the equation to calculate the tangent plane equation
tan_plane_eqn = z - p(3) == fxv*(x - p(1)) + fyv*(y - p(2));

% output
disp('equation of the tangent plane:'); disp(tan_plane_eqn);


%% Problem 14
clear; clc
% Consider the function f(x, y) = (x + 7*y)*sin(x*y).
% (a) Find the linearization L(x,y) of f at (1, π/2).
% (b) Using the linearization L(x,y) in part (a) above, find an 
% approximation for f(1.01, π/2 - 0.02).

syms x y;
f = (x + 7*y)*sin(x*y);
p = [1 pi/2];
delta_p = [1.01 pi/2 - 0.02];

% LINEARIZATION EQUATION:
% L(x, y) = f(a, b) + fx(a, b)*(x - a) + fy(a, b)*(y - b)
% (a)
% get partial derivatives
fx = diff(f, x);
fy = diff(f, y);
% get fx(a, b), fy(a, b), and f(a, b)
fxv = subs(fx, [x y], p);
fyv = subs(fy, [x y], p);
fv = subs(f, [x y], p);
% use equation
L = fv + fxv*(x - p(1)) + fyv*(y - p(2));
% (b)
% plug in delta p
Lv = subs(L, [x y], delta_p);
% output
disp('L(x, y):'); disp(L);
disp('L(1.01, π/2 - 0.02):'); disp(Lv);

%% Problem 15
clear; clc
% Given z = e^(x^2 - 3*y), find the differential dz.

syms x y;
z = exp(x^2 - 3*y);

% find partial derivatives
dzdx = diff(z, x);
dzdy = diff(z, y);
% combine into a single derivative
syms dx dy;
dz = dzdx*dx + dzdy*dy;
% output
fprintf('dz = %s dx + %s dy\n', dzdx, dzdy);

%% Problem 16
clear; clc
% Let w = 3*x^2*y*z + y*z^2
% (a) Find the differential dw
% (b) If x changes from 0.01, y changes from 1 to 0.98 and z changes from 1
% to 1.05, find dw. 

syms x y z dx dy dz;
w = 3*x^2*y*z + y*z^2;
p = [0 1 1];
dv = [0.01, 0.98, 1.05] - p;

% find partial derivatives
dwdx = diff(w, x);
dwdy = diff(w, y);
dwdz = diff(w, z);
% combine into a single derivative
dw = dwdx*dx + dwdy*dy + dwdz*dz;
% plug in the x, y, z values at point p and the changes dx, dy, dz.
dwv = subs(dw, [x y z dx dy dz], [p dv]);
% output
fprintf('dw = %s dx + %s dy + %s dz\n', dwdx, dwdy, dwdz);
fprintf('dwv = %s\n', dwv);

%% Problem 17
clear; clc
% Given z = x^2 - 3*x*y where x = t^2 - 1 and y = t^3 + t, find dz/dt at
% t = 1

syms x y t;
z = x^2 - 3*x*y;
tv = 1;

% CHAIN RULE FORMULA:
% dz/dt = (dz/dx)*(dx/dt) + (dz/dy)*(dy/dt)
% get partial derivatives of z
dzdx = diff(z, x);
dzdy = diff(z, y);
% define x and y
x = t^2 - 1;
y = t^3 + t;
% get derivatives of x and y
dxdt = diff(x, t);
dydt = diff(y, t);
% plug into chain rule formula
dzdt = dzdx*dxdt + dzdy*dydt;
% substitute in equations for x and y
dzdtv = subs(subs(dzdt), t, tv);
% output
fprintf('dz/dt|t = 1 = %s\n', dzdtv);

%% Problem 18
clear; clc
% Given w = x*y^2*z^3 where x = r*cos(s), y = r*sin(s) and z = r^2*s, find 
% ∂w/∂s at (r, s) = (1, π/2).

syms x y z r s
w = x*y^2*z^3;
rsv = [1 pi/2];

% CHAIN RULE FORMULA:
% ∂w/∂s = (∂w/∂x)*(∂x/∂s) + (∂w/∂y)*(∂y/∂s) + (∂w/∂z)*(∂z/∂s)
% get partial derivatives of w with respect to x, y, and z
dwdx = diff(w, x);
dwdy = diff(w, y);
dwdz = diff(w, z);
% define x and y
x = r*cos(s);
y = r*sin(s);
z = r^2*s;
% get partial derivatives of x, y, and z with respect to s
dxds = diff(x, s);
dyds = diff(y, s);
dzds = diff(z, s);
% apply chain rule formula to get ∂w/∂s
dwds = dwdx*dxds + dwdy*dyds + dwdz*dzds;
% substitute in the equations for x, y, and z
dwdss = subs(dwds);
% substitute the point provided in to get the final value
dwdssv = subs(dwdss, [r s], rsv);
% output
fprintf('∂w/∂s|(r, s) = (1, π/2) = %s\n', dwdssv);

%% Problem 19
clear; clc
% Each dimension of a rectangular box without lid is increasing at a 
% constant rate of 2 cm/s. Find the rate of change of the surface area of 
% the box when the dimensions of the box are 10x5x4, where 4 is the height 
% of the box.

syms l w h;
S = l*w + 2*l*h + 2*w*h;    % the box has no lid
dimensions = [10 5 4];
dldt = 2;
dwdt = 2;
dhdt = 2;

% CHAIN RULE FORMULA:
% dS/dt = (dSdl)*(dldt) + (dSdw)*(dwdt) + (dSdh)*(dhdt)
% get the partial derivatives of S with respect to l, w, and h
dSdl = diff(S, l);
dSdw = diff(S, w);
dSdh = diff(S, h);
% apply the chain rule formula
dSdt = dSdl*dldt + dSdw*dwdt + dSdh*dhdt;
% plug in the current dimensions of the box
dSdtv = subs(dSdt, [l w h], dimensions);
% output
fprintf('dS/dt: %s\n', dSdtv);

%% Problem 20
clear; clc
% Suppose the function z is defined implicitly by the equation
% x^2*y*z^2 + 5*z*y = 1 - x. Find ∂z/∂x.

syms x y z;
% rearrange F to have all terms on one side
F = x^2*y*z^2 + 5*z*y - 1 + x;

% IMPLICIT DIFFERENTIATION FORMULA:
% ∂z/∂x = -(∂F/∂x)/(∂F/∂z)
% get partial derivatives of F with respect to x and z
dFdx = diff(F, x);
dFdz = diff(F, z);
% apply formula
dzdx = -dFdx/dFdz;
% output
fprintf('∂z/∂x = %s\n', dzdx);




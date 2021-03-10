%% Problem 1
% Find an equation of the tangent plane to the given surface at the 
% specified point. z = 8*x^2 + y^2 − 9y,    (1, 4, −12)
clear; clc

syms x y z;

% declare plane equation and the specified point
p_eqn = z == 8*x^2 + y^2 - 9*y;
p = [1 4 -12];

% EQUATION OF TANGENT PLANE:
% z − z0 = fx(x0, y0)*(x − x0) + fy(x0, y0)*(y − y0)
% to satisfy this, we need to take the partial derivatives with respect to
% x and y. Where z = f(x, y)
f = rhs(p_eqn);
fx = diff(f, x);
fy = diff(f, y);
% get values by plugging in x and y from the specified point
fxv = double(subs(fx, [x y], [p(1) p(2)]));
fyv = double(subs(fy, [x y], [p(1) p(2)]));
% use the equation to calculate the tangent plane equation
tan_plane_eqn = z - p(3) == fxv*(x - p(1)) + fyv*(y - p(2));

% output
disp('equation of the tangent plane:'); disp(tan_plane_eqn);

%% Problem 2
% Find an equation of the tangent plane to the given surface at the 
% specified point. z = ln(x − 4*y),    (5, 1, 0)
clear; clc

syms x y z;

% declare plane equation and the specified point
p_eqn = z == log(x - 4*y);
p = [5 1 0];

% EQUATION OF TANGENT PLANE:
% z − z0 = fx(x0, y0)*(x − x0) + fy(x0, y0)*(y − y0)
% to satisfy this, we need to take the partial derivatives with respect to
% x and y. Where z = f(x, y)
f = rhs(p_eqn);
fx = diff(f, x);
fy = diff(f, y);
% get values by plugging in x and y from the specified point
fxv = double(subs(fx, [x y], [p(1) p(2)]));
fyv = double(subs(fy, [x y], [p(1) p(2)]));
% use the equation to calculate the tangent plane equation
tan_plane_eqn = z - p(3) == fxv*(x - p(1)) + fyv*(y - p(2));

% output
disp('equation of the tangent plane:'); disp(tan_plane_eqn);

%% Problem 3
% Explain why the function is differentiable at the given point.
% f(x, y) = 1 + x*ln(x*y − 5),    (3, 2)
% Find the linearization L(x, y) of f(x, y) at (3, 2).
clear; clc

syms x y;

% declare f(x, y) and the specified point
f = 1 + x*log(x*y - 5);
p = [3 2];
% take the partial derivatives with respect to x and y
fx = diff(f, x);
fy = diff(f, y);
% plug values in for x and y
fxv = double(subs(fx, [x y], p));
fyv = double(subs(fy, [x y], p));
% we can tell just by looking at 5 that it is only defined when the value
% inside the ln is greater than 5. That is, x*y > 5
% since both derivatives are continuous and x*y > 5, f is differentiable at
% the specified point
% LINEARIZATION EQUATION:
% L(x, y) = f(a, b) + fx(a, b)*(x - a) + fy(a, b)*(y - b)
% we need to calculate f(a, b)
fpv = double(subs(f, [x y], p));
% now we can calculate L(x, y)
L = fpv + fxv*(x - p(1)) + fyv*(y - p(2));

% output
disp('fx(x, y):'); disp(fx);
disp('fy(x, y):'); disp(fy);
disp('fx(p):'); disp(fxv);
disp('fy(p):'); disp(fyv);
disp('L(x, y):'); disp(L);

%% Problem 4
% Find the linear approximation of the function f(x, y, z) = 
% (x^2 + y^2 + z^2)^(1/2) at (2, 9, 6) and use it to approximate the number 
% (2.01^2 + 8.97^2 + 5.98^2)^(1/2). 
% (Round your answer to five decimal places.)
clear; clc

syms x y z;

% declare f(x, y, z) and the specified point
f = (x^2 + y^2 + z^2)^(1/2);
p = [2 9 6];
% declare the point we will approximate
ap = [2.01 8.97 5.98];
% we will need the partial derivatives with respect to each variable
fx = diff(f, x);
fy = diff(f, y);
fz = diff(f, z);
% we will also need the values of these derivatives at the specified point
fxv = subs(fx, [x y z], p);
fyv = subs(fy, [x y z], p);
fzv = subs(fz, [x y z], p);
% finally, we need the value of the function at the specified point
fpv = subs(f, [x y z], p);
% the linearization equation extends into n dimensions
% find L(x, y, z)
L = fpv + fxv*(x - p(1)) + fyv*(y - p(2)) + fzv*(z - p(3));
% plug in the approximated point to find a value
Lv = double(subs(L, [x y z], ap));

% output (the question asks to be rounded to 5 decimal places)
format long
disp('L(ap):'); disp(round(Lv, 5));
format short

%% Problem 5
% Find the differential of the function.
% z = e^(−9*x)*cos(5*(pi)*t)
clear; clc

syms x t z;

% declare z(x, y)
z = exp(-9*x)*cos(5*pi*t);
% take the partial derivatives with respect to x and t
zx = diff(z, x);
zt = diff(z, t);
% MULTIVARIATE DERIVATIVE FORMULA:
% (d/dz) = (d/dx)*dx + (d/xt)*dt
dz = zx + zt;

% output
disp('dz:'); disp(dz);

%% Problem 6
% Find the differential of the function.
% m = p^6*q^7
clear; clc

syms m p q;

% declare m(p, q)
m = p^6*q^7;
% take partial derivatives with respect to p and q
mp = diff(m, p);
mq = diff(m, q);
% MULTIVARIATE DERIVATIVE FORMULA:
% (d/dm) = (d/dp)*dp + (d/xq)*dq
dm = mp + mq;

% output
disp('dm:'); disp(dm);

%% Problem 7
% If z = 3*x^2 + y^2 and (x, y) changes from (1, 3) to (1.05, 3.1),
% compare the values of Δz and dz. 
% (Round your answers to four decimal places.)
clear; clc

syms z x y;

% declare z(x, y) and the specified points
z = 3*x^2 + y^2;
p = [1 3];
cp = [1.05 3.1];
% Δz is just the difference between the values of z at both points
zp = double(subs(z, [x y], p));
zcp = double(subs(z, [x y], cp));
change_in_z = abs(zcp - zp);
% SMALL CHANGE FORMULA:
% dz = zx(1, 2)*dx + zy(1, 2)*dy 
% we need to calculate dx and dy as the change in x and y coordinates
dx = cp(1) - p(1);
dy = cp(2) - p(2);
% we also need both partial derivatives at point p
zx = diff(z, x);
zy = diff(z, y);
zxp = double(subs(zx, [x y], p));
zyp = double(subs(zy, [x y], p));
% plug use the formula to calculate dz
dz = zxp*dx + zyp*dy;

% output
disp('Δz:'); disp(change_in_z);
disp('dz:'); disp(dz);

%% Problem 8
% The length and width of a rectangle are measured as 48 cm and 26 cm, 
% respectively, with an error in measurement of at most 0.1 cm in each. 
% Use differentials to estimate the maximum error in the calculated area 
% of the rectangle.
clear; clc

syms x y;

% declare a function a(x, y) for the area of the rectangle
a = x*y;
x_and_y = [48 26];
% get partial derivatives
ax = diff(a, x);
ay = diff(a, y);
% the error is at most 0.1 for l and w
dx = 0.1;
dy = 0.1;
% SMALL CHANGE FORMULA:
% da = axv*dx + ayv*dy
% need ax and ay at the length and width
axv = double(subs(ax, [x y], x_and_y));
ayv = double(subs(ay, [x y], x_and_y));
% plug values into the small change formula
da = axv*dx + ayv*dy;

% output
disp('da:'); disp(da);

%% Problem 9
% A model for the surface area of a human body is given by 
% S = 0.1095*w^0.425*h^0.725, where w is the weight (in pounds), h is the 
% height (in inches), and S is measured in square feet. If the errors in 
% measurement of w and h are at most 2%, use differentials to estimate the 
% maximum percentage error in the calculated surface area. 
% (Round your answer to one decimal place.)

% this one resists matlab calculations simplifying it.
% dw/w = 0.02, dh/h = 0.02.
% dS/S = ... = 0.425(dw/w) + 0.725(dh/h)
dS = 0.425*0.02 + 0.725*0.02;

% output
disp('dS:'); disp(dS);

%% Problem 10
% Suppose you need to know an equation of the tangent plane to a surface S 
% at the point P(3, 1, 4).
% You don't have an equation for S but you know that the curves
% r1(t)	 = 	[3 + 3*t, 1 − t^2, 4 − 5*t + t^2]
% r2(u)	 = 	[2 + u^2, 2*u^3 − 1, 2*u + 2]
% both lie on S. Find an equation of the tangent plane at P.
clear; clc

syms x y z t u;

% declare r1(t), r2(u), and point p
r1 = [3 + 3*t 1 - t^2 4 - 5*t + t^2];
r2 = [2 + u^2 2*u^3 - 1 2*u + 2];
p = [3 1 4];
% we need to find t such that r1(t) = p and u such that r2(u) = p
tv = double(solve([r1(1) == p(1) r1(2) == p(2) r1(3) == p(3)], t));
uv = double(solve([r2(1) == p(1) r2(2) == p(2) r2(3) == p(3)], u));

% r1'(t) and r2'(u) represent the tangent vectors of r1(t) and r2(u)
% their cross product creates a vector tangent to the plane
% we can use this vector and the point p to generate a plane equation

% take derivatives of r1 and r2
r1p = diff(r1, t);
r2p = diff(r2, u);
% plug the values determined for t and u in to find tangent vectors at p
r1pv = double(subs(r1p, t, tv));
r2pv = double(subs(r2p, u, uv));
% take their cross product to find the normal vector to the tangent plane
nv = cross(r1pv, r2pv);
% plug nv and p into plane equation to get result
% PLANE EQUATION:
% nv(1)*(x - p(1)) + nv(2)*(y - p(2)) nv(3)*(z - p(3)) == 0
p_eqn = nv(1)*(x - p(1)) + nv(2)*(y - p(2)) + nv(3)*(z - p(3)) == 0;

% output
disp('plane equation:'); disp(p_eqn);

%% Problem 1
% Use the Chain Rule to find dz/dt.
% z = x*y^9 − x^2*y,    x = t^2 + 1,    y = t^2 − 1
clear; clc

syms x y z t;

% declare z(x, y)
z = x*y^9 - x^2*y;
% CHAIN RULE FORMULA:
% dz/dt = (dz/dx)*(dx/dt) + (dz/dy)*(dy/dt)
% we need partial derivatives
dzdx = diff(z, x);
dzdy = diff(z, y);
% now we can define x(t) and y(t)
x = t^2 + 1;
y = t^2 - 1;
% we also need derivatives of x and y with respect to t
dxdt = diff(x, t);
dydt = diff(y, t);
% apply the formula
dzdt = dzdx*dxdt + dzdy*dydt;

% output
disp('dz/dt:'); disp(dzdt);

%% Problem 2
% Use the Chain Rule to find dw/dt.
% w = x*e^(y/z),    x = t^3,    y = 3 − t,    z = 7 + 2*t
clear; clc

syms w x y z t;

% declare w(x, y, z)
w = x*exp(y/z);
% CHAIN RULE FORMULA:
% dw/dt = (dw/dx)*(dx/dt) + (dw/dy)*(dy/dt) + (dw/dz)*(dz/dt)
% we need partial derivatives
dwdx = diff(w, x);
dwdy = diff(w, y);
dwdz = diff(w, z);
% now we can define x(t), y(t), z(t)
x = t^3;
y = 3 - t;
z = 7 + 2*t;
% we also need derivatives of x, y and z with respect to t
dxdt = diff(x, t);
dydt = diff(y, t);
dzdt = diff(z, t);
% apply the formula
dwdt = dwdx*dxdt + dwdy*dydt + dwdz*dzdt;

% output
disp('dw/dt:'); disp(dwdt);

%% Problem 3
% Use the Chain Rule to find ∂z/∂s and ∂z/∂t.
% z = (x − y)^7,    x = s^2*t,    y = s*t^2
clear; clc

syms x y z s t;

% declare z(x, y)
z = (x - y)^7;
% CHAIN RULE FORMULAE:
% ∂z/∂s = (∂z/∂x)*(∂x/∂s) + (∂z/∂y)*(∂y/∂s)
% ∂z/∂t = (∂z/∂x)*(∂x/∂t) + (∂z/∂y)*(∂y/∂t)
% we need partial derivatives
dzdx = diff(z, x);
dzdy = diff(z, y);
% now we can define x(t) and y(t)
x = s^2*t;
y = s*t^2;
% we need derivatives of x and y with respect to t and s
dxdt = diff(x, t);
dydt = diff(y, t);
dxds = diff(x, s);
dyds = diff(y, s);
% apply the formula
dzds = dzdx*dxds + dzdy*dyds;
dzdt = dzdx*dxdt + dzdy*dydt;

% output
disp('∂z/∂s:'); disp(dzds);
disp('∂z/∂t:'); disp(dzdt);

%% Problem 4
% Use the Chain Rule to find the indicated partial derivatives.
% w = x*y + y*z + z*x,    x = r*cos(Ø),    y = r*sin(Ø),    z = r*Ø;
% ∂w/∂r, ∂w/∂Ø when r = 4, Ø = π/2
clear; clc

syms w x y z r o;

% declare w(x, y, z) and the specified point
w = x*y + y*z + z*x;
p = [4, pi/2];
% CHAIN RULE FORMULAE:
% ∂w/∂r = (∂w/∂x)*(∂x/∂r) + (∂w/∂y)*(∂y/∂r) + (∂w/∂z)*(∂z/∂r)
% ∂w/∂Ø = (∂w/∂x)*(∂x/∂Ø) + (∂w/∂y)*(∂y/∂Ø) + (∂w/∂z)*(∂z/∂Ø)
% get partial derivatives of w(x, y, z)
dwdx = diff(w, x);
dwdy = diff(w, y);
dwdz = diff(w, z);
% declare x(r, Ø), y(r, Ø), and z(r, Ø)
x = r*cos(o);
y = r*sin(o);
z = r*o;
% get all 6 partial derivatives
dxdr = diff(x, r);
dxdo = diff(x, o);
dydr = diff(y, r);
dydo = diff(y, o);
dzdr = diff(z, r);
dzdo = diff(z, o);
% find ∂w/∂r and ∂w/∂Ø
dwdr = (dwdx)*(dxdr) + (dwdy)*(dydr) + (dwdz)*(dzdr);
dwdo = (dwdx)*(dxdo) + (dwdy)*(dydo) + (dwdz)*(dzdo);
% plug in the specified point
% find values of x, y, and z at r and Ø
x = subs(x, [r o], p);
y = subs(y, [r o], p);
z = subs(z, [r o], p);
% find values of ∂w/∂r and ∂w/∂o at x, y, z, r, and Ø
dwdrv = subs(subs(dwdr, [r o], p));
dwdov = subs(subs(dwdo, [r o], p));

% output
disp('(∂w/∂r)(4, π/2):'); disp(dwdrv);
disp('(∂w/∂o)(4, π/2):'); disp(dwdov);

%% Problem 5
% Use this equation to find dy/dx.
% 2*y*cos(x) = x^2 + y^2
clear; clc

syms x y;

% declare equation F(x, y) as equal to zero
F = 2*y*cos(x) - x^2 + - y^2;
% EQUATION:
% dy/dx = -(∂F/∂x)/(∂F/∂y) = -Fx/Fy
% we need partial derivatives with respect to x and y
dFdx = diff(F, x);
dFdy = diff(F, y);
% apply equation
dydx = -dFdx/dFdy;

% output
disp('dy/dx:'); disp(dydx);

%% Problem 6
% Use the equations to find ∂z/∂x and ∂z/∂y.
% e^z = 9*x*y*z
clear; clc

syms x y z;

% declare F(x, y, z) as equal to zero
F = exp(z) - 9*x*y*z;
% EQUATIONS:
% ∂z/∂x = -(∂F/∂x)/(∂F/∂z)
% ∂z/dy = -(∂F/∂y)/(∂F/∂z)
% we need partial derivatives with respect to x, y, and z
dFdx = diff(F, x);
dFdy = diff(F, y);
dFdz = diff(F, z);
% apply equations
dzdx = -dFdx/dFdz;
dzdy = -dFdy/dFdz;

% output
disp('∂z/∂x:'); disp(dzdx);
disp('∂z/∂y:'); disp(dzdy);

%% Problem 7
% The radius of a right circular cone is increasing at a rate of 1.9 in/s 
% while its height is decreasing at a rate of 2.4 in/s. At what rate is 
% the volume of the cone changing when the radius is 135 in. and the 
% is 116 in.?
clear; clc

syms r h;

% declare volume equation V(r, h)
V = (1/3)*pi*r^2*h;
% to find the rate of change of volume, we need dV/dt. Use chain rule
% CHAIN RULE FORMULA:
% dV/dt = (∂V/∂r)*(dr/dt) + (∂V/dh)(dh/dt)
% we are given values for dr/dt and dh/dt
drdt = 1.9;
dhdt = -2.4;
% take partial derivatives
dVdr = diff(V, r);
dVdh = diff(V, h);
% we are also given values for r and h
r = 135;
h = 116;
% apply the chain rule formula to get an equation for dV/dt
dVdt = dVdr*drdt + dVdh*dhdt;
% plug values of r and h in to get final answer
dVdtv = subs(dVdt);

% output
disp('dV/dt:'); disp(dVdtv);

%% Problem 8
% The length ℓ, width w, and height h of a box change with time. 
% At a certain instant the dimensions are ℓ = 2 m and w = h = 4 m,
% and ℓ and w are increasing at a rate of 3 m/s while h is decreasing 
% at a rate of 7 m/s. At that instant find the rates at which the 
% following quantities are changing.
% (a) The volume.
% (b) The surface area.
% (c) The length of the diagonal. (Round your answer to two decimal
% places.)
clear; clc

syms l h w;

% declare functions for the volume, surface area, and length of the
% diagonal
V = l*h*w;
S = 2*(l*w + l*h + h*w);
L = sqrt(l^2 + h^2 + w^2);
% this is solved with chain rule, I hope the formula is obvious by now
% we are given values for dl/dt, dw/dt, and dh/dt
dldt = 3;
dwdt = 3;
dhdt = -7;
% take all 9 partial derivatives
dVdl = diff(V, l); dVdh = diff(V, h); dVdw = diff(V, w);
dSdl = diff(S, l); dSdh = diff(S, h); dSdw = diff(S, w);
dLdl = diff(L, l); dLdh = diff(L, h); dLdw = diff(L, w);
% use chain rule to calculate complete derivatives
dVdt = dVdl*dldt + dVdh*dhdt + dVdw*dwdt;
dSdt = dSdl*dldt + dSdh*dhdt + dSdw*dwdt;
dLdt = dLdl*dldt + dLdh*dhdt + dLdw*dwdt;
% declare given values for l, h, and w
l = 2;
h = 4;
w = 4;
% substitute values into derivatives
dVdtv = double(subs(dVdt));
dSdtv = double(subs(dSdt));
dLdtv = double(subs(dLdt));

% output
disp('(a) The volume:'); disp(dVdtv);
disp('(b) The surface area:'); disp(dSdtv);
disp('(c) The length of the diagonal:'); disp(round(dLdtv, 2));

%% Problem 9
% Find the directional derivative of f at the given point in the direction 
% indicated by the angle Ø.
% f(x, y) = x*y^3 − x^2,    (1, 4),    Ø = π/3
clear; clc

syms x y o;

% declare f(x, y)
f = x*y^3 - x^2;
% take partial derivatives
fx = diff(f, x);
fy = diff(f, y);
% If u is a unit vector in the direction of Ø, then use this equation
% DIRECTIONAL DERIVATIVE EQUATION:
% Duf(x, y) = fx(x, y)*cos(Ø) + fy(x, y)*sin(Ø)
Duf = fx*cos(o) + fy*sin(o);
% declare points and plug them in
x = 1;
y = 4;
o = pi/3;
Dufv = subs(Duf);

% output
fprintf('The value of Duf(x, y) at Ø is: %s\n', Dufv);

%% Problem 10
% Consider the following.
% f(x, y) = x/y,    P(4, 1),   u = <3/5, 4/5>
% (a) find the gradient of f.
% (b) evaluate ∆f at p
% (c) find the rate of change of f at p in the direction of u
clear; clc

syms x y;

% declare f(x, y, z), p, and u
f = x/y;
p = [4 1];
u = [3/5 4/5];
% (a)
% GRADIENT EQUATION:
% ∆f(x, y) = <fx(x, y), fy(x, y)>
Df = [diff(f, x) diff(f, y)];
% (b)
Dfpv = subs(Df, [x y], p);
% (c)
% DIRECTIONAL DERIVATIVE EQUATION:
% Duf(x, y) = ∆f(x, y) . u
Dufv = dot(Dfpv, u);

% output
disp('∆f(x, y, z):'); disp(Df);
disp('∆f(p):'); disp(Dfpv);
disp('Duf(p):'); disp(Dufv);


%% Problem 11
% Consider the following.
% f(x, y, z) = x^2*y*z − x*y*z^8,    P(6, −1, 1),   u = <0, 4/5, −3/5>
% (a) find the gradient of f.
% (b) evaluate ∆f at p
% (c) find the rate of change of f at p in the direction of u
clear; clc

syms x y z;

% declare f(x, y, z), p, and u
f = x^2*y*z - x*y*z^8;
p = [6 -1 1];
u = [0 4/5 -3/5];
% (a)
% GRADIENT EQUATION:
% ∆f(x, y, z) = <fx(x, y, z), fy(x, y, z), fz(x, y, z)>
Df = [diff(f, x) diff(f, y) diff(f, z)];
% (b)
Dfpv = subs(Df, [x y z], p);
% (c)
% DIRECTIONAL DERIVATIVE EQUATION:
% Duf(x, y, z) = ∆f(x, y, z) . u
Dufv = dot(Dfpv, u);

% output
disp('∆f(x, y, z):'); disp(Df);
disp('∆f(p):'); disp(Dfpv);
disp('Duf(p):'); disp(Dufv);

%% Problem 12
% Find the directional derivative of the function at the given point in 
% the direction of the vector v.
% g(s, t) = s*t^(1/2),    (2, 4),    v = 2i − j
clear; clc

syms s t;

% declare g(s, t), p, and v
g = s*sqrt(t);
p = [2 4];
v = [2 -1];
% since v is not a unit vector, find a vector u that is
u = v/norm(v);
% start by getting the gradient
% GRADIENT EQUATION:
% ∆g(s, t) = <gs(s, t), gt(s, t)>
Dg = [diff(g, s) diff(g, t)];
% get the value of the gradient at point p
Dgpv = subs(Dg, [s t], p);
% DIRECTIONAL DERIVATIVE EQUATION:
% Dug(s, t) = ∆g(s, t) . u
Dugv = dot(Dgpv, u);

% output
disp('Dug(p):'); disp(Dugv);

%% Problem 13
% Find the directional derivative of the function at the given point in 
% the direction of the vector v.
% f(x, y, z) = x^2*y + y^2*z,    (2, 7, 9),    v = <2, −1, 2>
clear; clc

syms x y z;

% declare f(x, y, z), p, and v
f = x^2*y + y^2*z;
p = [2 7 9];
v = [2 -1 2];
% since v is not a unit vector, find a vector u that is
u = v/norm(v);
% start by getting the gradient
% GRADIENT EQUATION:
% ∆f(x, y, z) = <fx(x, y, z), fy(x, y, z), fz(x, y, z)>
Df = [diff(f, x) diff(f, y) diff(f, z)];
% get the value of the gradient at point p
Dfpv = subs(Df, [x y z], p);
% DIRECTIONAL DERIVATIVE EQUATION:
% Duf(x, y, z) = ∆f(x, y, z) . u
Dufv = dot(Dfpv, u);

% output
disp('Duf(p):'); disp(Dufv);

%% Problem 14
% Find the maximum rate of change of f at the given point and the 
% direction in which it occurs. f(x, y) = 8*y*sqrt(x),  (16, 7)
clear; clc

syms x y;

% declare f(x, y) and p
f = 8*y*sqrt(x);
p = [16 7];
% MAXIMUM RATE OF CHANGE EQUATION:
% mrc = norm(∆f(x, y));
% MAXIMUM RATE OF CHANGE DIRECTION EQUATION:
% mrcd = ∆f(x, y)/mrc
% get the gradient
Df = [diff(f, x) diff(f, y)];
% get the direction vector from above equation
mrcd = subs(Df, [x y], p);
% get the maximum rate of change from the above equation
mrc = norm(mrcd);
% get the unit direction vector
umrcd = mrcd./mrc;

% output
disp('The maximum rate of change:'); disp(mrc);
disp('The unit direction vector:'); disp(umrcd);

%% Problem 15
% Find all points at which the direction of fastest change of the function 
% f(x, y) = x^2 + y^2 − 2*x − 4*y is <1, 1>
% (Enter your answer as an equation.)
clear; clc

syms x y;

% declare f(x, y) and v
f = x^2 + y^2 - 2*x - 4*y;
v = [1 1];
% the direction in which the maximum rate of change of f(x, y) occurs at a
% point (a, b) is given by ∆f(a, b)
Df = [diff(f, x) diff(f, y)];
% we need to find all points (x, y) for which ∆f(x, y) = k*v
% essentially, this boils down to k = ∆f(1), k = ∆f(2), the solutions are
% when these values are equal
eqn = Df(1) == Df(2);
% simplify by solveing for y
simp_eqn = y == solve(eqn, y);

% output
disp('values that satisfy requirements are given by:'); disp(simp_eqn);

%% Problem 16
% Suppose that over a certain region of space the electrical potential V 
% is given by the following equation.
% V(x, y, z) = 3*x^2 − 5*x*y + x*y*z
% (a) Find the rate of change of the potential at P(5, 4, 7) in 
% the direction of the vector v = <1 1 -1>.
% (b) In which direction does V change most rapidly at P?
% (c) What is the maximum rate of change at P?
clear; clc

syms x y z;

% V(x, y, z)
V = 3*x^2 - 5*x*y + x*y*z;
% (a)
% gradient function
DV = [diff(V, x) diff(V, y) diff(V, z)];
% p and v
p = [5 4 7];
v = [1 1 -1];
% v isn't a unit vector, make one
u = v./norm(v);
% get gradient at p
DVp = subs(DV, [x y z], p);
% get magnitude in direction of u
DVpu = dot(DVp, u);
% (b)
% that's just the value of DV at p
% (c)
% that's just the magnitude of DV at p
DVpm = norm(DVp);

% output
disp('∆V(x, y, z):'); disp(DV);
disp('∆V(p):'); disp(DVp);
disp('|∆V(p)|:'); disp(DVpm);sq

%% Problem 17
% Find equations of the tangent plane and the normal line to the given 
% surface at the specified point.
% x + y + z = 6*exp(x*y*z),    (0, 0, 6)
% (a) the tangent plane
% (b) the normal line
clear; clc

syms x y z t;

% declare F(x, y, z) as equal to 0, p
F = x + y + z - 6*exp(x*y*z);
p = [0 0 6];
% (a)
% this means that x + y + z = 6*exp(x*y*z) is the level surface 
% F(x, y, z) = 0.
% find the gradient
DF = [diff(F, x) diff(F, y) diff(F, z)];
% find the gradient at p
DFp = subs(DF, [x y z], p);
% apparently this value represents a normal vector to the tangent plane
% use this and p to construct a plane equation
p_eqn = DFp(1)*(x - p(1)) + DFp(2)*(y - p(2)) + DFp(3)*(z - p(3)) == 0;
% (b)
% the normal line heads in the same direction as the tangent plane
nl = p + DFp*t;

% output
disp('The equation for the tangent plane:'); disp(p_eqn);
disp('Parametric equations for the normal line:'); disp(nl);
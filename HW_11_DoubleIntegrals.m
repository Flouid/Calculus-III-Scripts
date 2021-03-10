% Louis Keith
% 2/28/21
% MA 231 Homework 11
%% Problem 1
clear; clc
% Evaluate the double integral by first identifying it as the volume of a 
% solid. ∫∫[R](12 − 6y)dA, R = [0, 1] × [0, 1]

syms y;
R = [0 1];
f = 12 - 6*y;

% evaluate the double integral
di_f = int(int(f, y, R(1), R(2)), y, R(1), R(2));
% output
fprintf('Value: %s\n', di_f);

%% Problem 2
clear; clc
% Calculate the iterated integral. 
% ∫[-8, 8]∫[0, π/2](y + y^2*cos(x))dx dy

syms x y;
f = y + y^2*cos(x);
r1 = [0 pi/2];
r2 = [-8 8];

% evaluate integral
di_f = int(int(f, x, r1(1), r1(2)), y, r2(1), r2(2));
% output
fprintf('Value: %s\n', di_f);

%% Problem 3
clear; clc
% ∫∫[R](x*sec(y)^2 dA,    R = {(x, y) | 0 ≤ x ≤ 8, 0 ≤ y ≤ π/4}

syms x y;
f = x*sec(y)^2;
R = [0 8; 0 pi/4];

% evaluate integral
di_f = int(int(f, x, R(1, 1), R(1, 2)), y, R(2, 1), R(2, 2));
% output
fprintf('Value: %s\n', di_f);

%% Problem 4
clear; clc
% Find the volume of the solid in the first octant bounded by the parabolic
% cylinder z = 16 - x^2 and the plane y = 1.

syms x y;
z = 16 - x^2;

% The cylinder intersects the xy-plane along the line x = 4, so in the
% first octant, the solid lies below the surface z = 16 − x^2 and above the 
% rectangle R = [0, 4] × [0, 1] in the xy-plane.
R = [0 4; 0 1];

% take the double integral ∫∫[R](16 - x^2) dA,  R = [0, 1] x [0, 4]
volume = int(int(z, x, R(1, 1), R(1, 2)), y, R(2, 1), R(2, 2));
% output
fprintf('Volume: %s\n', volume);

%% Problem 5
clear; clc
% Evaluate the iterated integral.
% ∫[0 2]∫[0 y^2](x^2*y)dx dy

syms x y;
f = x^2*y;
R = [0 y^2; 0 2];

value = int(int(f, x, R(1, 1), R(1, 2)), y, R(2, 1), R(2, 2));
% output
fprintf('Value: %s\n', value);

%% Problem 6
clear; clc
% Consider the following.
% ∫∫[D](x) dA,    D is enclosed by the lines y = x, y = 0, x = 2

syms x y;
f = x;
Rt1 = [0 x; 0 2];
Rt2 = [y 2; 0 2];

% As a type I region, D lies between the lower boundary y = 0 and the upper 
% boundary y = x for 0 ≤ x ≤ 2, so D = {(x, y) | 0 ≤ x ≤ 2, 0 ≤ y ≤ x}.
% If we describe D as a type II region, D lies between the left boundary 
% x = y and the right boundary x = 2 for 0 ≤ y ≤ 2, so 
% D = {(x, y) | 0 ≤ y ≤ 2, y ≤ x ≤ 2}.
t1v = int(int(f, y, Rt1(1, 1), Rt1(1, 2)), x, Rt1(2, 1), Rt1(2, 2));
t2v = int(int(f, x, Rt2(1, 1), Rt2(1, 2)), y, Rt2(2, 1), Rt2(2, 2));
% output
fprintf('Type 1 Value: %s\n', t1v);
fprintf('Type 2 Value: %s\n', t2v);


%% Problem 7
clear; clc
% Set up iterated integrals for both orders of integration. Then evaluate 
% the double integral using the easier order.
% ∫∫[D](y) dA,  D is bounded by y = x - 56; x = y^2

syms x y;
f = y;

% The curves y = x - 6 or x = y + 6 and x = y^2 intersect when y + 6 == y^2
x1 = y + 56;
x2 = y^2;
y_sol = solve(x1 == x2, y);

% as a type 1 region, the lower and upper boundaries are described by
% differing curves at different values of x.
% as a type 2 region, the left and right regions are always descrbied by
% the same curves for differing values of y.

% we describe D as a type II region, D is enclosed by the left boundary
% x = y^2 and right boundary x = y + 6 for y_sol(1) ≤ y ≤ y_sol(2), so 
% D = {(x, y) | y_sol(1) ≤ y ≤ y_sol(2), y^2 ≤ x ≤ y + 6}
% ∫∫[D](y) dA,  R = [y^2, y + 6] x [y_sol(1), y_sol(2)];
R = [y^2 y + 56; y_sol(1) y_sol(2)];
area = int(int(f, x, R(1, 1), R(1, 2)), y, R(2, 1), R(2, 2));
% output
fprintf('area: %s\n', area);

%% Problem 8
clear; clc
% Evaluate the double integral.
% ∫∫[D](2*y^2)dA,   D is the triangular region with vertices 
% (0, 1), (1, 2), (4, 1)


syms x y;
f = 2*y^2;

% declare points
p1 = [0 1]; p2 = [1 2]; p3 = [4 1];
% since p1 and p3 are on the same y, they form the "base" of the triangle/
% create equations for lines from each side of the base to the top point
% p2.
m1 = (p2(2) - p1(2))/(p2(1) - p1(1));   % rise over run
l1 = y - p1(2) == m1*(x - p1(1));       % point-slope form
x1 = solve(l1, x);                      % solve for x
m2 = (p2(2) - p3(2))/(p2(1) - p3(1));
l2 = y - p3(2) == m2*(x - p3(1));
x2 = solve(l2, x);

% a type II region makes the most sense to use
% ∫[p2(1), p2(2)]∫[x1, x2](3*y^2)dx dy
area = int(int(f, x, x1, x2), y, p2(1), p2(2));
% output
fprintf('area: %s\n', area);

%% Problem 9
clear; clc
% Find the volume of the given solid. 
% Under the plane 7*x + 2*y - z = 0 and above the region enclosed by the
% parabolas y = x^2 and x = y^2.

syms x y;
z = 7*x + 2*y;  % plane equation
l1 = y == x^2;  % line 1
l2 = x == y^2;  % line 2

% solve both lines for y, only take positive solutions
y1 = solve(l1, y);
y2 = solve(l2, y); y2 = y2(2); % the positive one
% this can be thought of as a type I region where √x lies above x^2
[xp, yp] = solve(l1, l2, x <= 1, y <= 1);
% ∫[xp(1), xp(2)]∫[y1 y2](7*x + 2*y)dy dx
volume = int(int(z, y, y1, y2), x, xp(1), xp(2));
% output
fprintf('volume: %s\n', volume);

%% Problem 10
clear; clc
% Find the volume of the given solid. 
% Enclosed by the paraboloid z = x^2 + y^2 + 1 and the planes x = 0, y = 0,
% z = 0, and x + y = 5.

syms x y;
z = x^2 + y^2 + 1;
plane = x + y == 5;

% the boundaries x = 0, y = 0, and z = 0 are essentially saying this is
% only inside the first octant. 
l1 = solve(plane, y);
% get points
[x1, y1] = solve(x == 0, plane);
[x2, y2] = solve(y == 0, plane);
% evaluate as a type 1 region
% ∫[x1 x2]∫[y2 l1](x^2 + y^2 + 1)dy dx
volume = int(int(z, y, y2, l1), x, x1, x2);
% output
fprintf('volume: %s\n', volume);

%% Problem 11
clear; clc
% visual order of integration problem

%% Problem 12
clear; clc
% visual order of integration problem

%% Problem 13
clear; clc
% Evaluate the integral by reversing the order of integration.
% ∫[0 3]∫[3*y 9](11e^(x^2))dx dy

syms x y;
f = 11*exp(x^2);

% ∫[0 9]∫[0 x/3](11*exp(x^2))dy dx
volume = int(int(f, y, 0, x/3), x, 0, 9);
% output
fprintf('volume: %s\n', volume);

%% Problem 14
clear; clc
% Evaluate the integral by reversing the order of integration
% ∫[0 64]∫[y^(1/3) 4](8e^(x^4))dx dy

syms x y;
f = 8*exp(x^4);

% ∫[0 4]∫[0 x^3](8*exp(x^4))dy dx
volume = int(int(f, y, 0, x^3), x, 0, 4);
% output
fprintf('volume: %s\n', volume);
%% Problem 1
% Consider the point P(-1, 3, 7). Find the distance from P to the 
% (a) xz-plane
% (b) x-axis
clear; clc

p = [1 -3 7];

% the distance to the xz-plane is just absolute value of the y-coordinate
a = abs(p(2));
disp(a);

% the distance to the x-axis is just the square root of the yz coordinates
b = sqrt(p(2)^2 + p(3)^2);
disp('sqrt:'); disp(b^2);

%% Problem 2
% Find the distance between two points P(-2, 5, 4) and Q(1, 3, 2).
clear; clc

p = [-2 5 4];
q = [1 3 2];

% find the distance between each coordinate
pq = q - p;
% square root of the sum of the squares is the distance formula
d = norm(pq);
disp('sqrt:'); disp(d^2);

%% Problem 3
% Find the vector PQ where P and Q are the points (4, 9, 7) and (5, 3, 8)
% respectively
clear; clc

p = [4 9 7];
q = [5 3 8];

% calculate the vector
pq = q - p;
disp(pq);

%% Problem 4
% Let a = <6, -1, 4> and b = <2, 1, 0>. Find 2a - b.
clear; clc

a = [6 -1 4];
b = [2 1 0];

% calculate result
result = 2*a - b;
disp(result);

%% Problem 5
% Let a = <1, 2, -1>.
% (a) Find |a| (the magnitude of a).
% (b) Find a unit vector in the direction of a.
% (c) Find a vector with a magnitude of 3 in the opposite direction of a.
clear; clc

a = [1 2 -1];

% (a)
magnitude = norm(a);
disp('sqrt:'); disp(magnitude^2);

% (b)
% this is just a./|a|
% [1/sqrt(6), 2/sqrt(6), -1/sqrt(6)]

% (c)
% this is just -3*(a./|a|)
% [-3/sqrt(6) -6/sqrt(6) 3/sqrt(6)]

%% Problem 6
% Let a = <7, -2, 3> and b = <1, 1, 2>.
% (a) Find the dot product a.b.
% (b) Find the angle theta between a and b where 0 <= theta <= pi.
% (c) Are the two vectors orthogonal? "Y" or "N".
clear; clc

a = [7 -2 3];
b = [1 1 2];

% (a)
dp = dot(a, b);
disp(dp);

% (b)
% the dot product has two formulae: a(1)*b(1) + a(2)*b(2) + a(3)*b(3) = dp,
% and |a|*|b|*cos(theta) = dp. Rearrange to get dp/(|a|*|b|) = cos(theta)
disp(norm(a)^2);
disp(norm(b)^2);
% sqrt(6) and sqrt(62)
% cos_theta = 11/(sqrt(6)*sqrt(62))
% cos^(-1)(11/(sqrt(6)*sqrt(62))) = theta

% (c) 
% since the dot product does not equal zero, the answer is N

%% Problem 7
% A vector a of magnitude 2 makes angles alpha and beta with the positive
% x-axis and the positive y-axis respectively, where cos(alpha) = 0.4 and
% cos(beta) = 0.5.
% (a) Find the cosine of the angle gamme that a makes with the positive
% x-axis.
% (b) Find the vector a.
clear; clc

cos_alpha = 0.4;
cos_beta = 0.5;
magnitude_a = 2;

% (a) 
% cos_alpha is defined as a_1/|a|. This means a_1 = cos_alpha*|a|
a_1 = cos_alpha*magnitude_a;
% similarly, cos_beta is defined as a_2/|a|. Same logic applies
a_2 = cos_beta*magnitude_a;
% also, sqrt(a_1^2 + a_2^2 + a_3^2) = |a|. Rearrange and get
a_3 = sqrt(magnitude_a^2 - a_1^2 - a_2^2);
disp('sqrt:'); disp(a_3^2);
% finally, cos_gamma = a_3/|a|
% matlab won't output that cleanly, the answer is sqrt(59/25)/2

% (b)
a = [a_1 a_2 a_3];
disp(a);

%% Problem 8
% In the figure below, the coordinatew of points A and B are (6, 3, 5) and
% (1, 2, 3) respectively.
% (a) Find the vector OC.
% (b) Find the vector AC.
clear; clc;

a = [6 3 5];
b = [1 2 3];
o = [0 0 0];

% in the figure, it looks like C is the projection of A onto B
% note that proj_b_onto_a = (a.b/|a|)(a/|a|). This means that 
% proj_a_onto_b = (a.b/|b|)(b/|b|);
magnitude_b = norm(b);
a_dot_b = dot(a, b);
% we are only concerned in the squared magnitude of b
magnitude_b_squared = magnitude_b^2;
% proj_a_onto_b = (27/14)b
c = (a_dot_b/magnitude_b_squared)*b;
disp(c);
% c = [27/14 27/7 81/14];

% (a)
oc = c - 0;
disp(oc);

% (b)
ac = c - a;
disp(ac);

%% Problem 9
% Let a = <5, 1, -1> and b = <3, 0, -2>.
% (a) Find the cross product a x b.
% (b) Find the area of the triangle whose two sides are along the vectors a
% and b.
clear; clc

a = [5 1 -1];
b = [3 0 -2];

% (a)
a_cross_b = cross(a, b);
disp(a_cross_b);

% (b)
% the magnitude of the cross product defines the area of a parallelogram
% formed by two vectors. The triangle would simply be half of this value.
numerator = norm(a_cross_b);
area = numerator/2;
disp('sqrt:'); disp(numerator^2);
% answer is sqrt(62)/2

%% Problem 10
% Let a = 9i + j - 2k and b = 3i + j - k. The volume of the parallelopiped
% whose three adjacent edges are along the vectors a, b, and xk is 60
% units, where x is a positive real number. Find the value of x.
clear; clc

syms x;

a = [9 1 -2];
b = [3 1 -1];
c = [0 0 x];        % the vector xk
volume = 60;

% the volume of a parallelopiped is defined as the scalar triple product of
% the three vectors that define it. The scalar triple product is 
% (a x b).c
% calculate the scalar triple product
stp = dot(cross(a, b), c);
% set up an equation to make it equal the volume and solve for x
solution = solve(stp == volume, x);
disp('x = '); disp(solution);

%% Problem 11
% The line l passes through the point A(1, 2, 0) and is perpindicular to
% the plane p1 whose equation is given by 3x - y + 4z + 1 = 0.
% (a) Find the equation for l in parametric form. Use t for the parameter.
% (b) Find the linear equation of the plane p2 that passes through the
% point A and is parallel to the plane p1.
clear; clc

syms x y z t;

a = [1 2 0];
nv = [3 -1 4];      % the normal vector to the plane

% (a) 
% if the line is perpindicular to the plane, that means that we can think
% of it as beginning at point a and travelling in the diretion normal to
% the plane.
l = a + nv*t;
disp(l);

% (b)
% we can generate a plane equation using a normal vector and a point. Since
% p1 and p2 are parallel, they have the same normal vector. For a point p
% and normal vector nv, the equation is given by 
% nv(1)(x - p(1)) + nv(2)(y - p(2)) + nv(3)(z - p(3)) == 0
plane_eqn = nv(1)*(x - a(1)) + nv(2)*(y - a(2)) + nv(3)*(z - a(3)) == 0;
disp(plane_eqn);

%% Problem 12
% Find the linear equation of the plane that contains the y-axis and passes
% through the point (-1, 1, 3)
clear; clc

syms x y z;

p = [-1 1 3];

% if the plane contains the y-axis, then it must also contain the origin
o = [0 0 0];
% this means that the line through the origin and p lie in the plane
op = p - o;
% the y-axis is defined as a vector in the [0 1 0] direction, and also lies
% in the plane. Taking the cross product of these two lines should provide
% a normal vector to the plane
nv = cross(op, [0 1 0]);
% since the origin lies in the plane, we can use it and the normal vector
% to define its equation
plane_eqn = nv(1)*x + nv(2)*y + nv(3)*z == 0;
disp(plane_eqn);

%% Problem 13
% Find the shortest distance between the two lines given below
% (x - 3)/2 = y/1 = (z+2)/-2
% x/3 = (y-1)/2 = (z-2)/2
clear; clc

syms t;

% lets define these parametrically
l1 = [2*t + 7 t -2*t - 2];      % l1 = p + u*t
l2 = [3*t 2*t + 1 2*t + 2];     % l2 = q + v*t
% from these parametric equations, we can define as a point plus a
% direction vector times t. 
p = [7 0 -2];
u = [2 1 -2];
q = [0 1 2];
v = [3 2 2];
% apparently, the distance between two lines is defined as
% |pq.(u x v)|/|(u x v)| where u and v are the direction vectors and pq is
% the line generated by the points
pq = q - p;
u_cross_v = cross(u, v);
denominator = norm(u_cross_v)^2;
numerator = abs(dot(pq, u_cross_v));
disp('numerator:'); disp(numerator);
disp('denominator: sqrt'); disp(denominator);

%% Problem 14
% Find the coordinates of the point of intersection of the plane
% 2x + z + 2 = 0 and the line r = [1 0 -3] + t[2 1 1].
clear; clc

syms x y z t;

% declare the line r parametrically
r = [1 + 2*t t -3 + t];
% declare the plane equation
plane_eqn = 2*x + z + 2 == 0;
% substitue the parametrix equations in for x, y, and z.
parametric_plane_eqn = subs(plane_eqn, [x y z], r);
% solve for t
t = solve(parametric_plane_eqn, t);
% substitute t back into the line equation for a result
result = double(subs(r));
disp(result);

%% Problem 15
% Find the angle between the plane r.[5 1 -5] = 5 and the yz-plane.
clear; clc

syms x y z;

% so apparently that is a valid plane equation if r = [x y z];
p1_eqn = 5*x + y -5*z - 5 == 0;
p2_eqn = x == 0;
% the normal vectors for these
nv1 = [5 1 -5];
nv2 = [1 0 0];
% calculate the dot product
dp = dot(nv1, nv2);
% calculate the magnitudes of the normal vectors
magnitude_1_squared = norm(nv1)^2;
magnitude_2_squared = norm(nv2)^2;
disp(magnitude_1_squared);
disp(magnitude_2_squared);
% calculate the cos of theta
cos_theta = dp/(sqrt(magnitude_1_squared)*sqrt(magnitude_2_squared));
% the answer is 5/sqrt(51)

%% Problem 16
% The angle between a plane and a line is defined as the angle between the
% line and its projection on the plane. This, in turn, is equal to the
% complementary acute angle formed between the normal vector of the plane
% and the direction vector of the line. Find the angle between the plane 
% 7x + 2y + z + 6 = 0 and the line x/3 = y - 1 = (z - 2)/-1.
clear; clc

syms t;

% declare the normal vector for the plane
nv = [7 2 1];
% the line is given symmetrically, convert to parametric form
l = [3*t t + 1 -t + 2];
% declare the direction vector for the line
lv = [3 1 -1];
% find the dot product of the vectors
dp = dot(nv, lv);
% calculate the magnitudes
magnitude_1_squared = norm(nv)^2;
magnitude_2_squared = norm(lv)^2;
disp(magnitude_1_squared);
disp(magnitude_2_squared);
% calculate the cos of theta
cos_theta = dp/(sqrt(magnitude_1_squared)*sqrt(magnitude_2_squared));

%% Problems 17-24
% Look on page 837 in the textbook. There is a pattern. Memorize.
% x^2 + y^2 + z^2 = 1       ellipsoid
% z = x^2 + y^2             elliptic paraboloid
% z = x^2 - y^2             hyperbolic paraboloid
% z^2 = x^2 + y^2           cone
% x^2 + y^2 - z^2 = 1       hyperboloid of one sheet
% -x^2 -y^2 + z^2 = 1       hyperboloid of two sheets

%% Problem 25
% Choose a vector function that represents the curve of intersection of the
% cylinder y^2 + z^2 = 16 and the surface x + 3z = y^2.
clear; clc

syms x y z t;

circle_eqn = y^2 + z^2 == 16;
surface_eqn = x + 3*z == y^2;

% first, we parametrize the circle
% y = radius * sin(t)
% z = radius * cos(t)
radius = sqrt(double(rhs(circle_eqn)));
y = radius*sin(t);
z = radius*cos(t);
% substitute the parametric equations for y and z into the surface equation
parametric_surface_eqn = subs(surface_eqn);
% solve the parametric surface equation for x
x = solve(parametric_surface_eqn, x);
% use x, y, and z to make a vector function for the curve of intersection
vector = [x y z];
disp(vector);

%% Problem 26
% Find the derivative r'(t) where r(t) = [t^3 ln(7 + t) tan(t)].
clear; clc

syms t;

r = [t^3 log(7 + t) tan(t)];

% take the derivative
r_prime = diff(r, t);
disp(r_prime);

%% Problem 27
% Find the a unit tangent vector T(t) to the curve C given by 
% r(t) = [2 - t 6 + t^2 sin(6t)] at t = 0.
clear; clc

syms t;

T = [2 - t 6 + t^2 sin(6*t)];
% take the derivative
T_prime = diff(T, t);
% plug in t = 0
vector = double(subs(T_prime, t, 0));
disp(vector);
magnitude_squared = norm(vector)^2;
disp('sqrt:'); disp(magnitude_squared);



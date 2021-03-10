%% Problem 1
% Consider the point P(9, -7, 9). Find the distance from P to the 
% (a) xz-plane
% (b) x-axis
clear; clc

p = [9 -7 9];

% the distance to the xz-plane is just absolute value of the y-coordinate
a = abs(p(2));
% the distance to the x-axis is just the square root of the yz-coordinates
b = sqrt(p(2)^2 + p(3)^2);

% output
disp('(a):'); disp(a);
disp('(b):'); disp(subs(b));

%% Problem 2
% Find the distance between two points P(-2, 9, 4) and Q(1, 3, 2).
clear; clc

p = [-2 9 4];
q = [1 3 2];

% find the distance between each coordinate
pq = q - p;
% square root of the sum of the squares is the distance formula
d = norm(pq);

% output
disp(subs(d));

%% Problem 3
% Find the vector PQ where P and Q are the points (8, 5, 8) and (5, 3, 8)
% respectively
clear; clc

p = [8 5 8];
q = [5 3 8];

% calculate the vector
pq = q - p;

% output
disp(pq);

%% Problem 4
% Let a = <1, -1, 4> and b = <2, 1, 0>. Find 2a - b.
clear; clc

a = [1 -1 4];
b = [2 1 0];

% calculate result
result = 2*a - b;

% output
disp(result);

%% Problem 5
% Let a = <6, 2, -1>.
% (a) Find |a| (the magnitude of a).
% (b) Find a unit vector in the direction of a.
% (c) Find a vector with a magnitude of 3 in the opposite direction of a.
clear; clc

a = [6 2 -1];

% (a)
magnitude = norm(a);
% (b)
unit_vector = a./magnitude;
% (c)
vector = -3.*unit_vector;

% output
disp('(a):'); disp(subs(magnitude));
disp('(b):'); disp(subs(unit_vector));
disp('(c):'); disp(subs(vector));

%% Problem 6
% Let a = <2, -2, 3> and b = <1, 1, 2>.
% (a) Find the dot product a.b.
% (b) Find the angle theta between a and b where 0 <= theta <= pi.
% (c) Are the two vectors orthogonal? "Y" or "N".
clear; clc

a = [2 -2 3];
b = [1 1 2];

% (a)
dp = dot(a, b);
% (b)
% the dot product has two formulae: a(1)*b(1) + a(2)*b(2) + a(3)*b(3) = dp,
% and |a|*|b|*cos(theta) = dp. Rearrange to get dp/(|a|*|b|) = cos(theta)
cos_theta = dp/(norm(a)*norm(b));
% (c) 
% since the dot product does not equal zero, the answer is N

% output
disp('(a):'); disp(dp);
disp('(b):'); disp(subs(cos_theta));
disp('(c):'); disp('The dot product does not equal 0, thus the vectors are not orthogonal');

%% Problem 7
% A vector a of magnitude 2 makes angles alpha and beta with the positive
% x-axis and the positive y-axis respectively, where cos(alpha) = 0.8 and
% cos(beta) = 0.5.
% (a) Find the cosine of the angle gamme that a makes with the positive
% x-axis.
% (b) Find the vector a.
clear; clc

cos_alpha = 0.8;
cos_beta = 0.5;
magnitude_a = 2;

% (a) 
% cos_alpha is defined as a_1/|a|. This means a_1 = cos_alpha*|a|
a_1 = cos_alpha*magnitude_a;
% similarly, cos_beta is defined as a_2/|a|. Same logic applies
a_2 = cos_beta*magnitude_a;
% also, sqrt(a_1^2 + a_2^2 + a_3^2) = |a|. Rearrange and get
a_3 = sqrt(magnitude_a^2 - a_1^2 - a_2^2);
% finally, cos_gamma = a_3/|a|
cos_gamma = a_3/magnitude_a;
% (b)
a = [a_1 a_2 a_3];

% output
disp('(a):'); disp(subs(cos_gamma));
disp('(b):'); disp(subs(a));

%% Problem 8
% In the figure below, the coordinatew of points A and B are (2, 3, 5) and
% (1, 2, 3) respectively.
% (a) Find the vector OC.
% (b) Find the vector AC.
clear; clc;

a = [2 3 5];
b = [1 2 3];
o = [0 0 0];

% in the figure, it looks like C is the projection of A onto B
% note that proj_b_onto_a = (a.b/|a|)(a/|a|). This means that 
% proj_a_onto_b = (a.b/|b|)(b/|b|);
magnitude_b = norm(b);
a_dot_b = dot(a, b);
c = (a_dot_b/magnitude_b).*(b./magnitude_b);
% (a)
oc = c - o;
% (b)
ac = c - a;

% output
disp('(a):'); disp(subs(oc));
disp('(b):'); disp(subs(ac));

%% Problem 9
% Let a = <1, 1, -1> and b = <3, 0, -2>.
% (a) Find the cross product a x b.
% (b) Find the area of the triangle whose two sides are along the vectors a
% and b.
clear; clc

a = [1 1 -1];
b = [3 0 -2];

% (a)
a_cross_b = cross(a, b);
% (b)
% the magnitude of the cross product defines the area of a parallelogram
% formed by two vectors. The triangle would simply be half of this value.
area = (1/2)*norm(a_cross_b);

% output
disp('(a):'); disp(a_cross_b);
disp('(b):'); disp(subs(area));

%% Problem 10
% Let a = 14i + j - 2k and b = 3i + j - k. The volume of the parallelopiped
% whose three adjacent edges are along the vectors a, b, and xk is 60
% units, where x is a positive real number. Find the value of x.
clear; clc

syms x;

a = [14 1 -2];
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

% output
disp('x = '); disp(solution);

%% Problem 11
% The line l passes through the point A(1, 2, 0) and is perpindicular to
% the plane p1 whose equation is given by 3x - y + 8z + 1 = 0.
% (a) Find the equation for l in parametric form. Use t for the parameter.
% (b) Find the linear equation of the plane p2 that passes through the
% point A and is parallel to the plane p1.
clear; clc

syms x y z t;

a = [1 2 0];
nv = [3 -1 8];      % the normal vector to the plane

% (a) 
% if the line is perpindicular to the plane, that means that we can think
% of it as beginning at point a and travelling in the diretion normal to
% the plane.
l = a + nv*t;
% (b)
% we can generate a plane equation using a normal vector and a point. Since
% p1 and p2 are parallel, they have the same normal vector. For a point p
% and normal vector nv, the equation is given by 
% nv(1)(x - p(1)) + nv(2)(y - p(2)) + nv(3)(z - p(3)) == 0
plane_eqn = nv(1)*(x - a(1)) + nv(2)*(y - a(2)) + nv(3)*(z - a(3)) == 0;

% output
disp('(a):'); disp(l);
disp('(b):'); disp(plane_eqn);

%% Problem 12
% Find the linear equation of the plane that contains the y-axis and passes
% through the point (-1, 1, 7)
clear; clc

syms x y z;

p = [-1 1 7];

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

syms x y z t;

l1_sym = [(x - 3)/2 y (z + 2)/-2];
l2_sym = [x/3 (y - 1)/2 (z - 2)/2];

% we can convert to parametric form
l1 = [solve(l1_sym(1) == t, x) solve(l1_sym(2) == t, y) solve(l1_sym(3) == t, z)];
l2 = [solve(l2_sym(1) == t, x) solve(l2_sym(2) == t, y) solve(l2_sym(3) == t, z)];
% we need to extract a point and direction vector from each of these lines
% the points can be extracted by plugging in t = 0
l1_point = double(subs(l1, t, 0));
l2_point = double(subs(l2, t, 0));
% the direction vectors can be extracted by plugging in t = 1 and
% subtracting the point values
l1_dir = double(subs(l1, t, 1)) - l1_point;
l2_dir = double(subs(l2, t, 1)) - l2_point;
% the distance between two lines is defined as
% |pq.(u x v)|/|(u x v)| where u and v are the direction vectors and pq is
% the line generated by the points
pq = l2_point - l1_point; 
cp = cross(l1_dir, l2_dir);
distance = abs(dot(pq, cp))/norm(cp);

% output
disp('distance:'); disp(subs(distance));

%% Problem 14
% Find the coordinates of the point of intersection of the plane
% 6x + z + 2 = 0 and the line r = [1 0 -3] + t[2 1 1].
clear; clc

syms x y z t;

r_point = [1 0 -3];
r_dir = [2 1 1];
plane_eqn = 6*x + z + 2 == 0;

% calculate r
r = r_point + r_dir*t;
% substitute the parametric equations for r into the plane equation
plane_eqn_param = subs(plane_eqn, [x y z], r);
% solve for t
t = double(solve(plane_eqn_param, t));
% plug that t-value back into the equation for r
intersection = double(subs(r));

% output
disp('intersection:'); disp(subs(intersection));

%% Problem 15
% Find the angle between the plane r.[9 1 -5] = 5 and the yz-plane.
clear; clc

syms x y z;

r = [x y z];
r_dir = [9 1 -5];
intercept = 5;

% create the plane equation (dot product function return conjugates)
p1_eqn = r(1)*r_dir(1) + r(2)*r_dir(2) + r(3)*r_dir(3) == intercept;
% the yz plane is defined only in terms of x, and passes through the origin
p2_eqn = x == 0;
% the angle between the planes can be calculated by using the dot product
% of their normal vectors
% extract the normal vectors from each plane equation
nv1 = fliplr(double(coeffs(lhs(p1_eqn), [x y z]))); % needs a flip
nv2 = [double(coeffs(lhs(p2_eqn))) 0 0]; % two of the coefficients are 0
% calculate the dot product
dp = dot(nv1, nv2);
% using both of the formulae for the dot product,
% cos_theta = dp/(|nv1|*|nv2|)
cos_theta = dp/(norm(nv1)*norm(nv2));

% output
disp('cos(theta):'); disp(subs(cos_theta));

%% Problem 16
% The angle between a plane and a line is defined as the angle between the
% line and its projection on the plane. This, in turn, is equal to the
% complementary acute angle formed between the normal vector of the plane
% and the direction vector of the line. Find the angle between the plane 
% 4x + 2y + z + 6 = 0 and the line x/3 = y - 1 = (z - 2)/-1.
clear; clc

syms x y z t;

plane_eqn = 4*x + 2*y + z + 6 == 0;
l_sym = [x/3 y - 1 (z - 2)/-1];
intercept = 6;

% get the normal vector for the plane
% first, move the intercept to the right
plane_eqn = plane_eqn - intercept;
nv = fliplr(double(coeffs(lhs(plane_eqn))));
% convert the line to parametric form
l = [solve(l_sym(1) == t, x) solve(l_sym(2) == t, y) solve(l_sym(3) == t, z)];
% get the direction vector for the line
l_point = double(subs(l, t, 0));
l_dir = double(subs(l, t, 1)) - l_point;
% the same dot product and cos trick can be used to find the angle
dp = dot(nv, l_dir);
cos_theta = dp/(norm(nv)*norm(l_dir));

% output
disp('cos(theta):'); disp(subs(cos_theta));

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
% cylinder y^2 + z^2 = 49 and the surface x + 8z = y^2.
clear; clc

syms x y z t;

circle_eqn = y^2 + z^2 == 49;
surface_eqn = x + 8*z == y^2;

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

% output
disp('vector:'); disp(vector);

%% Problem 26
% Find the derivative r'(t) where r(t) = [t^3 ln(7 + t) tan(t)].
clear; clc

syms t;

r = [t^3 log(4 + t) tan(t)];

% take the derivative
r_prime = diff(r, t);
disp(r_prime);

%% Problem 27
% Find the a unit tangent vector T(t) to the curve C given by 
% r(t) = [2 - t 3 + t^2 sin(3t)] at t = 0.
clear; clc

syms t;

T = [2 - t 3 + t^2 sin(3*t)];
% take the derivative
T_prime = diff(T, t);
% plug in t = 0
vector = double(subs(T_prime, t, 0));
unit_vector = vector./norm(vector);

% output
disp('unit vector:'); disp(subs(unit_vector));

%% Various Exam Formulae

% VECTOR BETWEEN TWO POINTS FORMULA:
% <x2 - x1, y2 - y1, z2 - z1>

% VECTOR ADDITION
% r1 + r2 = <x1 + x2, y1 + y2, z1 + z2>

% SCALAR MULTIPLICATION
% a*r1 = <a*x1, a*y1, a*z1>

% UNIT VECTOR
% <x, y, z>/sqrt(x^2 + y^2 + z^2)

% DOT PRODUCT
% u . v = x1*x2 + y1*y2 + z1*z2
% u . v = |u|*|u|*cos(Ø)
% positive values mean 0 < Ø < 90
% negative values mean 90 < Ø < 180
% 0 means Ø = 90

% CROSS PRODUCT
% this would be a bitch to write down in matlab, its matrix determinants


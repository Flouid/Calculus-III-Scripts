%% Problem 2
clear
clc

syms t;

p = [8 -9 6];
v = [1 5 -2/3];

r = p + v*t;
disp(r);

%% Problem 3
clear
clc

syms t; 

p = [0, 14, -6];
v = [3 -2 8];

r = p + v*t;
disp(r);

%% Problem 4
clear
clc

syms t;

p = [1 0 4];
v = [1 2 1];

r = p + v*t;
disp(r);

%% Problem 5
clear
clc

syms t;

a = [0 1/2 1];
b = [5 1 -3];
ab = b - a;

r = a + ab*t;
disp(r);

%% Problem 6
clear
clc

syms t;

a = [2 -2 8];
v = [-1 2 -3];

r = a + v*t;
disp(r);

%% Problem 7
clear
clc

syms t s;

l1 = [(6 + 4*t) (8 - 2*t) (2 + 6*t)];
l2 = [(1 + 4*s) (3 - 2*s) (4 + 5*s)];

disp(find_intersection(l1, l2));

%% Problem 9
% Find an equation of the plane. The plane through the point (3, 7, 8) 
% and with normal vector 4i + j − k
clear
clc

syms x y z;

p = [3 7 8];
% normal vector
nv = [4 1 -1];

disp(get_plane_eqn(p, nv));

%% Problem 10
% Find an equation of the plane. The plane through the point (3, 0, 1)
% and perpendicular to the line x = 4t, y = 4 − t, z = 1 + 7t
clear
clc

syms x y z;

p = [3 0 1];
% normal vector
nv = [4 -1 7];

disp(get_plane_eqn(p, nv));

%% Problem 11
% Find an equation of the plane. The plane through the point  (8, −1, −5) 
% and parallel to the plane 6x − y − z = 6
clear
clc

syms x y z;

p = [8 -1 -5];
plane = 6*x - y - z;

% since the point is parallel to the plane, it has the same function but a
% different intercept. We can plug the new point in to calculate a new
% intercept
disp(subs(plane, {x y z}, p));

%% Problem 12
% Find an equation of the plane.
% The plane that passes through the line of intersection of the planes 
% x − z = 2 and y + 2z = 2 and is perpendicular to the plane x + y − 2z = 4
clear
clc

syms x y z;

n1 = [1 0 -1]; % x - z = 2
n2 = [0 1 2];  % y + 2z = 2

% set z = 0 and find that (2, 2, 0) is a point on the line of intersection
p = [2, 2, 0];

% the direction of the line of intersection is the cross product of the
% planes
v1 = cross_product_3d(n1, n2);
% the plane is defined as x + y -2z = 4
v2 = [1 1 -2];

% the normal vector to the intersection is the cross product
nv = cross_product_3d(v1, v2);

disp(get_plane_eqn(p, nv));

%% Problem 14 
% Find the point at which the line intersects the given plane.
clear
clc

syms x y z t;

% parametric equations for each direction
x = 3 - 3*t;
y = 4*t;
z = 1 + t;
% equation for the plane
plane = x + 2*y - z == 6;

% substitute parametric equations into the plane equation
parametric_plane = subs(plane);
% solve the substitued equation for t
t = solve(parametric_plane);
% substitute that t back into the parametric equations to get a solution
solution = subs([x y z]);
% display the solution
disp(subs(solution));

%% Problem 15
% Find the cosine of the angle between the planes x + y + z = 0 and 
% x + 4y + 3z = 2.
clear
clc

nv1 = [1 1 1]; % normal vector for plane 1
nv2 = [1 4 3]; % normal vector for plane 2

% the angle between the two planes is the same as the angle between their
% normal vectors. 
% recall a.b = |a||b|cos(theta). => cos(theta) = a.b/(|a||b|)
format rat;
cos_theta = dot(nv1, nv2)/(norm(nv1)*norm(nv2));
disp(cos_theta);

%% Problem 16
% Find symmetric equations for the line of intersection of the planes.
% z = 4x − y − 13, z = 6x + 3y − 15
clear
clc

% declare normal vectors
nv1 = [4 -1 -1]; % normal vector for plane 1
nv2 = [6 3 -1];  % normal vector for plane 2
% the direction of the intersecting line is perpindicular to both planes.
v = double(cross_product_3d(nv1, nv2));

syms x y z;

% declare plane equations
p1eqn = z == 4*x - y - 13;
p2eqn = z == 6*x + 3*y - 15;
% if we substitute z = 0 into both equations, it reduces to a solvable
% system of equations.
simplified_eqns = [subs(p1eqn, z, 0), subs(p2eqn, z, 0)];
% solve the system
solution = solve(simplified_eqns, [x y]);
% this point is on the line of intersection and crosses the z-axis.
p = [double(solution.x) double(solution.y) 0];

% the symmetric equations are given by (x - x0)/v0
sym_eqn1 = (x - p(1))/v(1);

% outputting the first should be sufficient to determine the right choice.
disp(sym_eqn1);

%% Problem 17
% Find an equation for the plane consisting of all points that are 
% equidistant from the points (3, 0, −2) and (5, 10, 0).
clear
clc

% declare points of interest
p1 = [3 0 -2];
p2 = [5 10 0];

syms x y z;

% the distance formula to a point x y z, calculated for p1 and p2.
d1 = sqrt((x - p1(1))^2 + (y - p1(2))^2 + (z - p1(3))^2);
d2 = sqrt((x - p2(1))^2 + (y - p2(2))^2 + (z - p2(3))^2);
% the plane is defined as any point where these distances are equal.
plane = d1 == d2;
% simplify and display the equation. It is solved for x.
result = x == solve(plane);
disp(result);

%% Problem 18
% Find an equation of the plane with x-intercept a, y-intercept b, 
% and z-intercept c.
clear
clc

syms a b c;

% declare points of interest
xint = [a 0 0];
yint = [0 b 0];
zint = [0 0 c];

% now we need two vectors in the plane for a cross product. We can get
% these by drawing lines from the x-int to the y-int and z-int.
xy = yint - xint;
xz = zint - xint;
% find the normal vector using the cross product.
nv = cross(xy, xz);
% we can use the x-int as the point, and then we get a plane equation
disp(get_plane_eqn(xint, nv));

%% Problem 19
% Find parametric equations for the line through the point (0, 1, 2)
% that is parallel to the plane x + y + z = 1 and perpendicular to the line 
% x = 1 + t, y = 1 − t, z = 2t. (Use the parameter t.)
clear
clc

syms t;

% declare the point of interest.
p = [0 1 2];
% the normal vector of the plane.
nv = [1 1 1];
% the other vector.
v = [1 -1 2];

% the direction vector for the required line is the cross product of the
% normal vector of the plane and the other vector.
dv = cross(nv, v);

% finally, parametrize it on t and add it to the starting point
solution = p + dv*t;
disp(solution);

%% Problem 20
% Let P be a point not on the line L that passes through the points Q and R. 
% The distance d from the point P to the line L is d = |a × b|/|a|. 
% Where a = QR and b = QP. Use the above formula to find the distance from 
% the point to the given line. (0, 2, 3);    x = 3t, y = 5 − 2t, z = 3 + t
clear
clc

syms x y z t; 

% declare parametric equations for the line.
x = 3*t;
y = 5 - 2*t;
z = 3 + t;
% declare the point p.
p = [0 2 3];

% to get a vector a = qr, we can plug t = 0 and t = 1 into the line
% equation and subtract them.
q = double([subs(x, t, 0) subs(y, t, 0) subs(z, t, 0)]);
r = double([subs(x, t, 1) subs(y, t, 1) subs(z, t, 1)]);
a = r - q;
% to get a vector b = qp, we can do the same.
b = p - q;

% now we can calculate the distance using the formula above.
d = norm(cross(a, b))/norm(a);
disp(d);

%% Problem 21
% Find the distance between the given parallel planes.
% 4z = 6y − 6x,    6z = 3 − 9x + 9y
clear
clc

syms x y z;

% declare plane equations
p1eqn = 4*z == 6*y - 6*x;
p2eqn = 6*z == 3 - 9*x + 9*y;

% we can demonstrate that the first plane passes through the origin by
% plugging x = y = 0 into its equation.
p = double([0, 0, subs(p1eqn, [x y], [0 0])]);

% now we can calculate the distance from p to the second plane.
nv2 = [9 9 -6]; % the normal vector for the second plane will be useful.
% plug the origin into the plane 2 equation.
const = double(subs(rhs(p2eqn), [x y], [0 0]));
d = const/norm(nv2);
disp(d);

%% Problem 22
% Let L1 be the line through the origin and the point (7, 0, −1).
% Let L2 be the line through the points (1, −1, 1) and (6, 1, 7).
% Find the distance between L1 and L2.
clear
clc

syms x y z t s;

% declare a point on each line
p1 = [0 0 0];
p2 = [1 -1 1];
% declare direction vectors for lines 1 and 2.
dv1 = [7 0 -1] - p1;
dv2 = [6 1 7] - p2;
% declare parametric equations for l1 and l2. 
param_l1 = [p1(1) + dv1(1)*t, p1(2) + dv1(2)*t, p1(3) + dv1(3)*t];
param_l2 = [p2(1) + dv2(1)*s, p2(2) + dv2(2)*s, p2(3) + dv2(3)*s];

% check if the direction vectors are the same. If so, the vectors are
% parallel.
disp(isequal(dv1, dv2));
% since they aren't, check if they intersect.
param_eqns = [param_l1; param_l2];
solution = solve(param_eqns, [t s]);
disp(islogical(solution.t));
% since they are neither parallel nor intersect, they must be skew. We can
% view the lines as existing in two parallel planes, whose normal vector is
% given as the cross product of the lines.
nv = cross(dv1, dv2);
% additionally, p1 lies on the first plane, and p2 on the second plane. We
% can use this to create a pair of plane equations. Both equal 0.
p1eqn = nv(1)*x + nv(2)*y + nv(3)*z;
% the first plane passes through the origin, so it has no intercept. Since
% the second plane is parallel, the only difference is it's intercept. 
intercept = nv(1)*p2(1) + nv(2)*p2(2) + nv(3)*p2(3);
p2eqn = p1eqn - intercept;
% finally, the distance is given by the absolute value of the intercept
% over the magnitude of the normal vector. 
d = intercept/norm(nv);
disp(d);
% since matlab won't spit this out as a radical, the numbeer under the
% square root is the sum of the squares of the normal vector.
disp(sum(nv.^2));
% final answer: 63/sqrt(2409)




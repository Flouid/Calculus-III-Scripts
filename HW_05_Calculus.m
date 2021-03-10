%% Problem 1
% Consider the given vector equation. r(t) = (2t − 2, t^2 + 3)
% (a) Find r'(t).
% (b) Sketch the plane curve together with the position vector r(t) and 
% the tangent vector r'(t) for the given value of t = 2.
clear; clc

syms t;

% define r
r = [2*t - 2 t^2 + 3];
% take the derivative
r_prime = diff(r, t);
disp(r_prime);

%% Problem 2
% Find the unit tangent vector T(t) at the point with the given value of 
% the parameter t. r(t) = cos(t)i + 8tj + 3 sin(2t)k,    t = 0
clear; clc

syms t;

% define r
r = [cos(t) 8*t 3*sin(2*t)];
% take the derivative
r_prime = diff(r, t);
disp(r_prime);
% plug in t = 0;
solution = double(subs(r_prime, t, 0));
disp(solution);
% they want a unit vector, so find the magnitude
magnitude_squared = norm(solution)^2;
disp('sqrt:'); disp(magnitude_squared);

%% Problem 3
% If r(t) = (3t, 3t^2, 2t^3), find r'(t), T(1), r''(t), and r'(t) × r ''(t)
clear; clc

syms t;

% define r
r = [3*t 3*t^2 2*t^3];
% take the first and second derivatives
r_prime = diff(r, t);
disp('r_prime:'); disp(r_prime);
r_double_prime = diff(r_prime, t);
disp('r_double_prime:'); disp(r_double_prime);
% plug in t = 1 to get T(1)
T_1 = double(subs(r_prime, t, 1));
disp('T(1):'); disp(T_1);
magnitude_squared = norm(T_1)^2;
disp('divided by sqrt:'); disp(magnitude_squared);
% find the cross product of the first and second derivatives
cp = cross(r_prime, r_double_prime);
disp('cross product:'); disp(cp);

%% Problem 4
% Find parametric equations for the tangent line to the curve with the 
% given parametric equations at the specified point.
% x = e^(−4t)*cos(4t), y = e^(−4t)*sin(4t), z = e^(−4t); (1, 0, 1)
clear; clc

syms t;

% declare the curve as a vector
r = [exp(-4*t)*cos(4*t) exp(-4*t)*sin(4*t) exp(-4*t)];
% declare the point
p = [1 0 1];
% find the value of t that makes r = p.
t_value = double(solve(r(1) == p(1), t));
% take the derivative of r
r_prime = diff(r, t);
% plug our value of t in to get the tangent vector
tangent_vector = subs(r_prime, t, t_value);
% parametric solution will start at the point and do in direction of the
% tangent vector
solution = p + tangent_vector*t;
disp('solution:'); disp(solution);

%% Problem 5
% At what point do the curves r1(t) = (t, 5 − t, 35 + t^2) and 
% r2(s) = (7 − s, s − 2, s^2) intersect?
% Find their angle of intersection, theta, correct to the nearest degree.
clear; clc

syms s t;

% declare r1 and r2
r1 = [t 5 - t 35 + t^2];
r2 = [7 - s s - 2 s^2];
% we need to find t and s such that r1(t) = r2(s)
param_eqns = r1 == r2;
% solve for t and s
solutions = solve(param_eqns);
% plug t back into r1 to get the intersection point
intersection_point = subs(r1, t, solutions.t);
disp('intersection point:'); disp(intersection_point);

% to find the angle of intersection, we need to find the angle between the
% tangent vectors at the intersection point. First, we need their
% derivatives
r1_prime = diff(r1);
r2_prime = diff(r2);
% plug our earlier solutions for t and s into r1_prime and r2_prime to get
% their tangent vectors at the intersection point
r1_tv = double(subs(r1_prime, t, solutions.t));
r2_tv = double(subs(r2_prime, s, solutions.s));
% from earlier homeworks, we know that the cos(theta) = u.v/(|u||v|)
cos_theta = dot(r1_tv, r2_tv)/(norm(r1_tv)*norm(r2_tv));
% take the inverse cosine in degrees to get the solution
theta = acosd(cos_theta);
disp('theta:'); disp(theta);

%% Problem 6
% Evaluate the integral. integral from 0 to 2 of (4t i − t^3 j + 4*t^7 k)
clear; clc

syms t;

% declare the vector
r = [4*t -t^3 4*t^7];
% take it's integral
R = int(r);
% evaluate from 0 to 2
solution = double(subs(R, t, 2) - subs(R, t, 0));
disp('solution:'); disp(solution);

%% Problem 7
% Find the length of the curve. r(t) = 7t, 3 cos(t), 3 sin(t),  −4 ≤ t ≤ 4

% this is done on the ipad, unless I wanna figure out trig identities in
% matlab its not worth doing here

%% Problem 8
% Consider the following vector function. r(t) = (5t, (1/2)t^2, t^2)
% (a) Find the unit tangent and unit normal vectors T(t) and N(t).
% (b) Use this formula to find the curvature. k(t) = |T'(t)|/|r'(t)|

% this was also done on the ipad, is hell

%% Problem 9
% Find the curvature of r(t) = (5t, t^2, t^3) at the point (5, 1, 1).
clear; clc

syms t;

% this can be done by using an equation k(t) = |r'(t) x r''(t)|/|r'(1)|^3

% define r and p
r = [5*t t^2 t^3];
p = [5 1 1];
% find t that satisfies r(t) = (5, 1, 1)
t = double(solve(r(1) == p(1)));
% find the first and second derivatives of r
r_p = diff(r);
r_pp = diff(r_p);
% plug the t value in to find values for r'(t) and r''(t)
r_pv = double(subs(r_p));
r_ppv = double(subs(r_pp));
% find the cross product
cp = cross(r_pv, r_ppv);
% find the magnitudes of the cross product and the first derivative value
cp_magnitude_squared = sum(cp.^2);
rpv_magnitude_squared = sum(r_pv.^2);
% output
disp('sqrt:'); disp(cp_magnitude_squared);
disp('divided by');
disp('sqrt:'); disp(rpv_magnitude_squared^3);

%% Problem 10
% Use this formula to find the curvature. y = 5x^4
% k(x) = |f''(x)|/|1 + f'(x)^2|^(3/2)
clear; clc

syms x;

% declare the equation
y = 5*x^4;
% find it's first and second derivatives
y_p = diff(y);
y_pp = diff(y_p);
% find k(x)
k = abs(y_pp)/abs(1 + y_p^2)^(3/2);
% output
disp('k(x):'); disp(k);

%% Problem 11
% Find vectors T, N, and B at the given point.
% r(t) = (t^2, (2/3)t^3, t),    (1, -2/3, -1)
clear; clc

syms t;

% declare r and p
r = [t^2 (2/3)*t^3 t];
p = [1 -2/3 -1];

% T = r'(t)/|r'(t)|
% find t such that r(t) = p
t = double(solve(r == p));
% find r'(t)
rp = diff(r);
% find the magnitude of r'(t) squared
rpm = sqrt(sum(rp.^2));
% create an expression for T
T = rp./rpm;
% plug in the value for t to get an answer
Tv = double(subs(T));
% output
disp('T:'); disp(Tv);

% N = T'(t)/|T'(t)|
% find T'(t)
Tp = diff(T);
% find |T'(t)|
Tpm = sqrt(sum(Tp.^2));
% create and expression for N
N = Tp./Tpm;
% plug a value in for t to get an answer
Nv = double(subs(N));
% output
disp('N:'); disp(Nv);

% B = T x N
B = cross(T, N);
% plug in a value for t to get an answer
Bv = double(subs(B));
% output
disp('B:'); disp(Bv);

%% Problem 12
% At what point on the curve x = t^3, y = 9t, z = t^4 is the normal plane 
% parallel to the plane 6x + 18y − 8z = 4?
clear; clc

syms t;

% this problem effectively asks where r'(t) is parallel to the normal
% vector to the provided plane

% define r and the normal vector
r = [t^3 9*t t^4];
nv = [6 18 -8];
% find r'(t)
rp = diff(r);
% find t such that r'(t) is parallel to nv
% this happens when t = -1
t = -1;
tpv = double(subs(rp));
% plug that t back into r to find the point
rt = double(subs(r));
disp('r(t):'); disp(rt);
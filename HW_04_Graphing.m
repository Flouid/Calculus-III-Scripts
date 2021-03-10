%% Problem 14
% Find an equation for the surface consisting of all points that are 
% equidistant from the point (−2, 0, 0) and the plane x = 2.
clear; clc

syms x y z;

% let the point p be an arbitary point (x, y, z)
p = [x y z];
% declare the starting point and the plane equation coefficients
q = [-2 0 0];
cf = [1 0 0 -2];        % coefficients vector
% we can caluclate the distance from p to q
distance_p_to_q = norm(q - p);
% recall the distance from a plane formula: 
% D = |ax + by + cz + d|/sqrt(a^2 + b^2 + c^2)
% since the coefficients for the plane equation are [1 0 0 -3],
% plug the coefficients into the plane distance equation
numerator = abs(cf(1)*x + cf(2)*y + cf(3)*z + cf(4));
denominator = sqrt(cf(1)^2 + cf(2)^2 + cf(3)^2);
distance_p_to_plane = numerator/denominator;
% the surface as defined as any point where these two are equal
surface_eqn = distance_p_to_q == distance_p_to_plane;
disp(surface_eqn);
% this simplifies down to -8x = y^2 + z^2

%% Problem 15
% Find the limit. 
% lim t→0 ( e^(−6t)i + t^2/sin(t)^2j + sin(4t)k )
clear; clc

syms t;

% well the first term goes to e^0 = 1

% the second term is 0/0, so L'Hopital it into 2t/(2*sin(x)*cos(x)
% use an identity to convert to 2t/sin(2x). Still 0/0
% now it is 2/2cos(2x), which can now substitute to 2/2 = 1

% the third term evaluates to 0

disp([1 1 0]);

%% Problem 16
% Sketch the curve with the given vector equation. Indicate with an arrow 
% the direction in which t increases. r(t) = (sin(4t), t).
clear; clc

syms y t;

% define r
r = [sin(4*t) t];
% we can substitute y into x
x = sin(4*y);

%% Problem 17
% Sketch the curve with the given vector equation. Indicate with an arrow 
% the direction in which t increases. r(t) = (t^2 − 8, t).
clear; clc

syms y t;

% define r
r = [t^2 - 8, t];
% we can substitute y into x
x = y^2 - 8;

%% Problem 20
% At what points does the curve r(t) = ti + (6t − t^2)k intersect the 
% paraboloid z = x^2 + y^2? (If an answer does not exist, enter DNE.)
clear; clc

syms x y z t;

% declare curve r and the paraboloid equation
r = [t 0 6*t - t^2];
paraboloid = z == x^2 + y^2;
% substitute parametric equations into the paraboloid
param_parab = subs(paraboloid, [x y z], r);
% solve for t
t_values = double(solve(param_parab, t));
% substitute the solutions into r for both solutions
disp(subs(r, t_values(1)));
disp(subs(r, t_values(2)));

%% Problem 21
% Two particles travel along the space curves r1(t) = (t, t^2, t^3)
% r2(t) = (1 + 2t, 1 + 6t, 1 + 14t). Find the points at which their 
% paths intersect. (If an answer does not exist, enter DNE.)
% Find the time(s) when the particles collide. 
% (Enter your answers as a comma-separated list. If an answer does not 
% exist, enter DNE.)
clear; clc

syms s t;

% define r1 and r2
r1 = [t t^2 t^3];
r2 = [1 + 2*t 1 + 6*t 1 + 14*t];
% see if a value of t satisfies both equations
t = double(solve(r1(1) == r2(1), t));
% check if it satisfies both equations
disp(isequal(subs(r1, t), subs(r2, t)));
% since it doesn't, they don't collide. 

% we need to find t and s such that r1(t) = r2(s)
r2 = subs(r2, s);
param_eqns = r1 == r2;
% solve for t and s
solutions = solve(param_eqns);
t_values = double(solutions.t);
s_values = double(solutions.s);
% plug both t's into r1 to find the points of intersection
disp(subs(r1, t_values(1)));
disp(subs(r1, t_values(2)));




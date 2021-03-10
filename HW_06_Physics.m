%% Problem 1
% Find the velocity, acceleration, and speed of a particle with the given 
% position function. r(t) = (-(1/2)t^2, t)
clear; clc

syms t;

% define r(t)
rt = [-0.5*t^2 t];
% v(t) is the derivative of r(t)
vt = diff(rt);
disp('v(t):'); disp(vt);
% a(t) is the derivateive of v(t)
at = diff(vt);
disp('a(t):'); disp(at);
% speed of the particle |v(t)| is the magnitude of the velocity vector
st = sqrt(sum(vt.^2));
disp('|v(t)|:'); disp(st);

%% Problem 2
% Find the velocity, acceleration, and speed of a particle with the given 
% position function. r(t)  = (e^t)i + (e^(4t))j
clear; clc

syms t;

% define r(t)
rt = [exp(t) exp(4*t)];
% v(t) is the derivative of r(t)
vt = diff(rt);
disp('v(t):'); disp(vt);
% a(t) is the derivateive of v(t)
at = diff(vt);
disp('a(t):'); disp(at);
% speed of the particle |v(t)| is the magnitude of the velocity vector
st = sqrt(sum(vt.^2));
disp('|v(t)|:'); disp(st);

%% Problem 3
% (a) Find the position vector of a particle that has the given acceleration 
% and the specified initial velocity and position.
% a(t) = (11t)i + (e^t)j + (e^(-t))k,    v(0) = k,    r(0) = j + k
% (b) On your own using a computer, graph the path of the particle.
clear; clc

syms t C;

% (a)
% define a(t), v(0), r(0)
at = [11*t exp(t) exp(-t)];
v0 = [0 0 1];
r0 = [0 1 1];
% v(t) is the integral of a(t) + C, where v(0) = k
vt = int(at);
% now that we have the general integral, we need to solve for C
vt0 = subs(vt, t, 0) + C == v0;
C = double([solve(vt0(1), C) solve(vt0(2), C) solve(vt0(3), C)]);
% now we can use C to create a specific v(t)
vt = vt + C;
% we can repeat this process for r(t)
clear C; syms C
% r(t) is the integral of v(t) + C, where r(0) = j + k
rt = int(vt);
% now that we have the general integral, we need to solve for C
rt0 = subs(rt, t, 0) + C == r0;
C = double([solve(rt0(1), C) solve(rt0(2), C) solve(rt0(3), C)]);
% now we can use C to create a specific r(t)
rt = rt + C;
disp('r(t):'); disp(rt);

% (b)
figure(1)
fplot3(rt(1), rt(2), rt(3));
xlabel('x');
ylabel('y');
zlabel('z');

%% Problem 4
% What force is required so that a particle of mass m has the position 
% function r(t) = (t^3)i + (2t^2)j + (t^3)k?
clear; clc

syms m t;

% define r(t)
rt = [t^3 2*t^2 t^3];
% F(t) = m*a(t). We need to find a(t). As before v(t) is the derivative of 
% r(t) and a(t) if the derivative of v(t)
vt = diff(rt);
at = diff(vt);
% calculate F(t)
Ft = at*m;
% output
disp('F(t):'); disp(Ft);

%% Problem 5
% Find the tangential and normal components of the acceleration vector.
% r(t) = cos(t)i + sin(t)j + (t)k
clear; clc

syms t;

% define r(t)
rt = [cos(t) sin(t) t];

% a(t)_T = (r'(t).r''(t))/|r'(t)|
% so we need to find r'(t) = v(t), r''(t) = a(t), and s(t) = |r'(t)|
vt = diff(rt);
at = diff(vt);
% we can calulate the numerator
numerator = double(vt(1)*at(1) + vt(2)*at(2) + vt(3)*at(3));
% and the denominator
denominator = sqrt(sum(vt.^2));
% we notice that the identity sin(x)^2 + cos(x)^2 = 1 is present
denominator = double(subs(denominator, cos(t)^2 + sin(t)^2, 1));
% output
disp('a(t)_T:');
disp(numerator);
disp('divided by'); 
disp('sqrt:'); disp(denominator^2);

% a(t)_N = |r'(t) x r''(t)|/|r'(t)|
% we need to find the magnitude of the cross product of v(t) and a(t)
cp = cross(vt, at);
% again, recognize the identity
cp = subs(cp, cos(t)^2 + sin(t)^2, 1);
% find it's magnitude
cpm = sqrt(sum(cp.^2));
% yet again, recognize the identity 
cpm = double(subs(cpm, cos(t)^2 + sin(t)^2, 1));
% output
disp('a(t)_N:');
disp(cpm/denominator);

%% Problem 6
% A manufacturer has modeled its yearly production function P (the monetary 
% value of its entire production in millions of dollars) as a Cobb-Douglas 
% function P(L, K) = 1.47*L^(0.65)*K^(0.35)
% where L is the number of labor hours (in thousands) and K is the invested 
% capital (in millions of dollars). Find P(110, 40) and interpret it. 
% (Round your answers to one decimal place.)
clear; clc

syms L K;

% define P(L, K)
P = (1.47*L^(0.65))*K^(0.35);
% plug in P(110, 40)
result = double(subs(P, [L K], [110 40]));
% output
disp('P(110, 40):'); disp(result);

%% Problem 7
% Let g(x, y) = cos(x + 3y).
% (a) Evaluate g(6, −2).
% (b) Find the domain of g.
% (c) Find the range of g.
clear; clc

syms x y; 

% define g(x, y)
g = cos(x + 3*y);
% plug in g(6, -2)
result = double(subs(g, [x y], [6 -2]));








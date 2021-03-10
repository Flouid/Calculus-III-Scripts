%% Problem 1
% Consider the space curve given by r(t) = <cos(2t), t, sin(2t)>
% Find the unit tangent vector T at the point (1, π, 0).
clear; clc

syms t

% declare values
r = [cos(2*t) t sin(2*t)];
p = [1 pi 0];
% find t that satisfies r(t) = p
tv = solve(r(2) == p(2), t);
% find r'(t)
rp = diff(r, t);
% calculate |r'(t)|
rpm = sqrt(rp(1)^2 + rp(2)^2 + rp(3)^2);
% make necessary simplifying substitution
rpms = subs(rpm, cos(2*t)^2 + sin(2*t)^2, 1);
% plug t into both
rpv = subs(rp, t, tv);
rpmv = subs(rpms, t, tv);

% output
disp(rpv/rpmv);
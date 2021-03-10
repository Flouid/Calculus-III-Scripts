%% Problem 1
% Find the local maximum and minimum values and saddle point(s) of the 
% function. If you have three-dimensional graphing software, graph 
% the function with a domain and viewpoint that reveal all the important 
% aspects of the function. (Enter your answers as a comma-separated list. 
% If an answer does not exist, enter DNE.)
% f(x, y) = x^2 + x*y + y^2 + 7*y
clear; clc

syms x y;

% f(x, y)
f = x^2 + x*y + y^2 + 7*y;
minmaxsaddle(f);

%% Problem 2
% Find the local maximum and minimum values and saddle point(s) of the 
% function. If you have three-dimensional graphing software, graph 
% the function with a domain and viewpoint that reveal all the important 
% aspects of the function. (Enter your answers as a comma-separated list. 
% If an answer does not exist, enter DNE.)
% f(x, y) = x*y − 4*x − 4*y − x^2 − y^2
clear; clc

syms x y; 

% f(x, y)
f = x*y - 4*x - 4*y - x^2 - y^2;
minmaxsaddle(f);

%% Problem 3
% Find the local maximum and minimum values and saddle point(s) of the 
% function. If you have three-dimensional graphing software, graph 
% the function with a domain and viewpoint that reveal all the important 
% aspects of the function. (Enter your answers as a comma-separated list. 
% If an answer does not exist, enter DNE.)
% f(x, y) = 3 - x^4 + 2*x^2 - y^2
clear; clc

syms x y;

% f(x, y)
f = 3 - x^4 + 2*x^2 - y^2;
minmaxsaddle(f);

%% Problem 4
% Find the local maximum and minimum values and saddle point(s) of the 
% function. If you have three-dimensional graphing software, graph 
% the function with a domain and viewpoint that reveal all the important 
% aspects of the function. (Enter your answers as a comma-separated list. 
% If an answer does not exist, enter DNE.)
% f(x, y) = x^3 - 3*x + 3*x*y^2
clear; clc

syms x y;

% f(x, y)
f = x^3 - 3*x + 3*x*y^2;
minmaxsaddle(f);

%% Problem 5
% Find the absolute maximum and minimum values of f on the set D.
% f(x, y) = x^2 + y^2 + x^2*y + 5,
% D = {(x, y) | |x| ≤ 1, |y| ≤ 1}
clear; clc

syms x y;

% f(x, y)
f = x^2 + y^2 + x^2*y + 5;
conditions = [abs(x) <= 1 abs(y) <= 1];
minmaxsaddle(f, conditions);
 
%% Problem 6
% Find the absolute maximum and minimum values of f on the set D.
% f(x, y) = x*y^2 + 9,    D = {(x, y) | x ≥ 0, y ≥ 0, x^2 + y^2 ≤ 3}
clear; clc

syms x y;

% f(x, y)
f = x*y^2 + 9;
conditions = [x >= 0 y >= 0 x^2 + y^2 <= 3];
minmaxsaddle(f, conditions);










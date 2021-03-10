%% Problem 15
clear
clc

a = [4, 5, 0];
b = [1, 0 7];

disp(cross_product_3d(a, b))

%% Problem 16
clear
clc

syms t

a = [t, 2, 1/t];
b = [t^2, t^2, 1];

disp(cross_product_3d(a, b))

%% Problem 17
clear
clc

a = [1, 0, -2];
b = [0, 1, 1];

cross_product = cross_product_3d(a, b);

disp(cross_product)

%% Problem 24
clear
clc

a = [2, -1, 4];
b = [5, 2, 1];

disp(cross_product_3d(a, b))
disp(cross_product_3d(b, a))

%% Problem 25
clear
clc

p = [0 -2 0];
q = [6 1 -3];
r = [5 3 1];

pq = q - p;
pr = r - p;
result = cross_product_3d(pq, pr);

disp(result);
disp(norm(result));

%% Problem 26
clear
clc

a = [1 2 4];
b = [-1 1 4];
c = [5 1 2];

disp(scalar_triple_product(a, b, c));

%% Problem 27
clear
clc

u = [1 5 -2];
v = [4 -1 0];
w = [6 9 -4];

disp(scalar_triple_product(u, v, w));

%% Problem 28
clear
clc

a = [1 2 3];
b = [2 -4 8];
c = [6 0 1];
d = [5 6 -4];

ab = b - a;
ac = c - a;
ad = d - a;

disp(scalar_triple_product(ab, ac, ad));
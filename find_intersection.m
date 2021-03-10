function intersection = find_intersection(v1, v2)
    syms t s;

    eqn1 = v1(1) == v2(1);
    eqn2 = v1(2) == v2(2);
    eqn3 = v1(3) == v2(3);

    [A, B] = equationsToMatrix([eqn1, eqn2, eqn3], [t s]);
    intersection = linsolve(A, B);
end
function cross_product = cross_product_3d(a, b)
    % perform calculation with symbolic variables and determinant math
    syms I J K;
    unit_vector = [I J K];
    symbolic_result = det([unit_vector; a; b]);
    
    % convert to a vector and return
    I = 1; J = 1; K = 1;
    cross_product = subs(subs(children(symbolic_result)));
end
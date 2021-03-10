function plane_eqn = get_plane_eqn(p, d)
    syms x y z;
    
    plane_eqn = d(1)*(x - p(1)) + d(2)*(y - p(2)) + d(3)*(z - p(3)) == 0;
end
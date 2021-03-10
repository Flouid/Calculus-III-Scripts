function result = scalar_triple_product(a, b, c)
    result = norm(dot(a, cross_product_3d(b, c)));
end
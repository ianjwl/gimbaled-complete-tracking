function cross_product_equivalent_matrix = hat_map(vector)
%HAT_MAP Computes the cross product equivalent matrix for a given vector.
[n,m] = size(vector);
if n~=3||m~=1
    error('Error in hat_map, line 5: Input vector was of incorrect dimensions.  It must be 3x1.')
end
cross_product_equivalent_matrix = ...
    [       0      -vector(3)   vector(2);
        vector(3)      0       -vector(1)
       -vector(2)   vector(1)       0       ];
end


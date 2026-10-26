function T = mrp_kinematic_matrix(mrp)

if any(isinf(mrp))
    T = inf(3,3);
else
    T = 1/4*((1-norm(mrp)^2)*eye(3) + 2*hat_map(mrp) + 2*(mrp*mrp'));
end

end
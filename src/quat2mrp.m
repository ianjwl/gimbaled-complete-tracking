function mrp = quat2mrp(q)

if size(q,2) == 4
    flip_order = true;
    q = q./vecnorm(q.').';
    q = q';
else 
    flip_order = false;
    q = q/norm(q);
end

mrp = q(2:4,:)./(1+q(1,:));

if flip_order
    mrp = mrp';
end

end
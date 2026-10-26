function value = theta_bar(alpha,rho)
% rho = U1/r

value = inf(size(alpha));

mask = alpha >= sqrt(1+rho^2);
value(mask) = asin(alpha(mask).^(-1));

mask = (alpha >= 1 - rho)&(alpha < sqrt(1+rho^2));
value(mask) = acos((alpha(mask).^2+rho.^2-1)./(2*rho*alpha(mask)));

end
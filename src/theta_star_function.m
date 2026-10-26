function [value,alpha_2pi,alpha_pi0,alpha_pi1,alpha_plus] = theta_star_function(alpha,rho)
if (rho >= 1)||(rho < 0); error('rho needs to be in [0,1)!'); end

base = (1-rho)/2;
alpha_2pi = base/4;
alpha_pi0 = base/2;
alpha_pi1 = 3/4*base;
alpha_plus = 5/4*base;

I = intssat(alpha, alpha_pi1, alpha_plus);
D = alpha - I;                      % s - lstarplus
D(alpha >= alpha_plus) = base;      % exact there; avoids cancellation for large alpha
vartheta = 2*atan2(sqrt(D.*(1 - D)), sqrt(I.*(1 + I)));
vartheta(alpha <= alpha_pi1) = pi;  % degenerate triangle, including alpha = 0

sat = ssat(alpha,alpha_2pi,alpha_pi0);

value = vartheta.*sat + 2*pi*(1-sat);

end
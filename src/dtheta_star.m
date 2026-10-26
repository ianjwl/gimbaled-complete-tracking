function dvalue = dtheta_star(alpha,rho)
%DTHETA_STAR  d/dalpha of theta_star(alpha,rho), elementwise in alpha.
%   DTHETA_STAR(alpha,rho) returns the first derivative of the target
%   attitude bound theta_star with respect to alpha = ||a*||/r. The inner
%   loop needs it to form ddt_theta_star = grad_theta_star*ddt_alpha, which
%   in turn feeds ddt_n_star and the ddt_omega_star terms.
%
%   Construction. theta_star is the outer zeta blend of the intersection
%   angle vartheta with the 2*pi plateau. Both pieces are differentiated
%   analytically:
%       I  = Z(alpha; alpha_pi1, alpha_plus)        inner ramp
%       D  = alpha - I, pinned to base for alpha >= alpha_plus
%       vartheta   = 2*atan2(sqrt(D.*(1-D)), sqrt(I.*(1+I)))
%       theta_star = vartheta*sat + 2*pi*(1-sat),  sat = zeta(alpha;alpha_2pi,alpha_pi0)
%   so that
%       theta_star' = vartheta'*sat + (vartheta - 2*pi)*dsat.
%
%   The derivative of vartheta is carried in cot(vartheta/2) form,
%       cot(vartheta/2) = sqrt((1+I)./D./(1-D)) =: g,
%       vartheta' = -2*g'/(1+g^2),
%   which avoids differencing the two nearly equal atan2 arguments near
%   alpha = 0. vartheta is identically pi wherever I = 0, so its derivative
%   vanishes there; the mask (alpha > alpha_pi1) & (I > 0) restricts the
%   closed form to the region where it is nonsingular, and the complementary
%   entries are left at zero, which is the continuous extension. On the mask,
%   D and Q = D.*(1-D) stay bounded away from zero, so no cancellation
%   arises. The mask is also what keeps dk = vm./(2*k) well defined: both
%   numerator and denominator tend to zero together as alpha -> alpha_pi1+,
%   and the ratio stays finite.
%
%   Domain. rho = U1/r is accepted on [0,1). rho = 1 would collapse the
%   transition interval [alpha_pi1, alpha_plus], and rho < 0 is not physical.
%   rho = 0 is the U_1 = 0 case; the breakpoints alpha_2pi = (1-rho)/8,
%   alpha_pi0 = (1-rho)/4, alpha_pi1 = 3(1-rho)/8, alpha_plus =
%   5(1-rho)/8 stay distinct there, so no degenerate (a,b) pair reaches this
%   function.
%
%   Inputs:
%       alpha   normalized actuation norm ||a*||/r, any array shape
%       rho     U1/r, scalar in [0,1)
%
%   Outputs:
%       dvalue  d/dalpha of theta_star(alpha,rho), same shape as alpha
%
%   Requires: intssat, ssat, dssat.
%
%   Called by INNER_LOOP_CONTROLLER, MULTIROTOR_JUMP_GUARD, and
%   MULTIROTOR_LYAPUNOV_FUNCTION.

if (rho >= 1)||(rho < 0); error('rho needs to be between 0 and 1!'); end
base = (1-rho)/2;
alpha_2pi = base/4;
alpha_pi0 = base/2;
alpha_pi1 = 3/4*base;
alpha_plus = 5/4*base;

% inner saturation on [alpha_pi1, alpha_plus]
I = intssat(alpha, alpha_pi1, alpha_plus);
v = ssat(alpha, alpha_pi1, alpha_plus);
D = alpha - I;                      % s - lstarplus
D(alpha >= alpha_plus) = base;      % exact there

vartheta = 2*atan2(sqrt(D.*(1 - D)), sqrt(I.*(1 + I)));
vartheta(alpha <= alpha_pi1) = pi;  % degenerate triangle, including alpha = 0

% vartheta == pi identically where I == 0, so its derivative is zero there
dvartheta = zeros(size(alpha));
m = (alpha > alpha_pi1) & (I > 0);
if any(m(:))
    Im = I(m); vm = v(m); Dm = D(m);
    k  = sqrt(Im);                      % sqrt(I)
    dk = vm./(2*k);                     % ratio of two small, accurately computed numbers
    Q  = Dm.*(1 - Dm);                  % >= alpha_pi1*(1 - alpha_pi1) > 0
    dQ = (1 - vm).*(1 - 2*Dm);
    L  = 0.5*(vm./(1 + Im) - dQ./Q);    % h'/h
    h  = sqrt((1 + Im)./Q);
    g  = k.*h;                          % cot(vartheta/2)
    dg = h.*(dk + k.*L);
    dvartheta(m) = -2*dg./(1 + g.^2);
end

% outer blend
sat  = ssat(alpha, alpha_2pi, alpha_pi0);
dsat = dssat(alpha, alpha_2pi, alpha_pi0);
dvalue = dvartheta.*sat + (vartheta - 2*pi).*dsat;
end
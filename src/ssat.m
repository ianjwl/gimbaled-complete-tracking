function y = ssat(x,a,b)
%SSAT  Smooth step function zeta(x;a,b), the building block of the zeta/Z pair.
%   y = SSAT(x,a,b) returns the C-infinity step of Appendix A,
%       zeta(x;a,b) = zeta_0((x-a)./(b-a)),
%       zeta_0(u)   = varpi(u)./(varpi(u) + varpi(1-u)),
%       varpi(u)    = 0 for u <= 0, exp(-1./u) for u > 0,
%   evaluated elementwise. zeta is 0 for x <= a, 1 for x >= b, and rises
%   smoothly, with all derivatives vanishing at both ends, across a < x < b.
%   The intermediate variable s is the normalized coordinate (x-a)./(b-a),
%   i.e. the paper's u; it shares a name with, but is not, the paper's raw
%   argument s.
%
%   Calling conventions, matched across the family (dssat, ddssat, dddssat,
%   ddddssat):
%       ssat(x)      ->  s = x,            transition on 0 < x < 1
%       ssat(x,a)    ->  s = x./a,         transition on 0 < x < a
%       ssat(x,a,b)  ->  s = (x-a)./(b-a)
%   The 3-argument form also accepts b < a, in which case zeta falls smoothly
%   from 1 to 0 rather than rising, exploiting the identity
%       1 - zeta(s;s0,s1) = zeta(s;s1,s0).
%
%   Degenerate interval: Appendix A leaves zeta undefined when s0 = s1
%   because the normalization 1./(b-a) is not defined. Here the convention
%       ssat(x,a,a) = sign(x-a)
%   is adopted, so zeta becomes the unit step about the single point a. This
%   is reachable at rho = 0, where alpha_0 = alpha_1 = 0; the feedback
%   schedule kappa_circ then uses zeta(alpha;alpha_0,alpha_1) directly as
%   grad_kappa_o, and uses its integral Z (see INTSSAT) for kappa_circ
%   itself.
%
%   Numerical notes: phi and psi are formed separately, so the numerator and
%   the denominator are each computed directly and neither is obtained as
%   1 - p; the ratio therefore stays accurate in both tails. A trailing guard
%   replaces any NaN in the result, choosing 1 or 0 by comparing the
%   normalized coordinate against the midpoint 0.5; this covers non-finite or
%   otherwise unusable input without perturbing the analytic branches. The
%   explicit assignments at x = -Inf and x = +Inf make behavior at these
%   points explicit.  
%
%   Inputs:
%       x      raw argument, any array shape
%       a, b   transition endpoints, scalars
%
%   Outputs:
%       y      zeta(x;a,b), same shape as x
%
%   Called by THETA_STAR_FUNCTION, DTHETA_STAR, DDTHETA_STAR,
%   INNER_LOOP_CONTROLLER, MULTIROTOR_JUMP_GUARD,
%   MULTIROTOR_LYAPUNOV_FUNCTION, and the local JSSAT helper of
%   ALTERNATE_REFERENCE_FUNCTION.

if nargin == 1
    s = x;
elseif nargin == 2
    s = x./a;
elseif nargin == 3
    s = (x-a)./(b-a);
    if a == b
        y = sign(x-b);
        return
    end
else
    error('Bad number of arguments in ssat().')
end

phi = exp(-1./s).*(s>0);
psi = exp(-1./(1-s)).*((1-s)>0);
y = phi./(phi + psi);

y(x==-Inf) = 0;
y(x==Inf) = 1;

if any(isnan(y))
    % disp('y is NaN!')
    y(isnan(y)&(s>0.5)) = 1;
    y(isnan(y)&(s<=0.5)) = 0;
end

end


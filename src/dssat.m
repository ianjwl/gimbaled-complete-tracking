function dy_by_dx = dssat(x,a,b)
%DSSAT  First derivative of ssat(x,a,b) with respect to x.
%   Companion to ssat.m, ddssat.m, dddssat.m, ddddssat.m. Same calling
%   conventions:
%       dssat(x)      ->  s = x,             ds/dx = 1
%       dssat(x,a)    ->  s = x./a,          ds/dx = 1/a
%       dssat(x,a,b)  ->  s = (x-a)./(b-a),  ds/dx = 1/(b-a)
%
%   IMPLEMENTATION NOTE
%   Because s(x) is affine, the chain rule collapses to a single factor,
%       dy/dx = (dy/ds)(ds/dx).
%   The closed form used for dy/ds is the derivative of phi./(phi+psi) with
%   phi = exp(-1/s) and psi = exp(-1/(1-s)), written as one rational
%   expression in those exponentials rather than assembled from the quotient
%   rule term by term. This gets the whole derivative in a single pass and
%   avoids cancellation between separately formed pieces. 
%
%   dy/ds is zeroed outside the support 0 < s < 1, where the transition has
%   either not begun or already completed, so dy/dx vanishes off (a,b) as it
%   must. NOTE: dssat and ddssat state their derivatives as single rational
%   closed forms in the exponentials, whereas dddssat and ddddssat build
%   theirs from the logistic polynomials L1..L4 in (p,q) = (y,1-y). The two
%   representations are equivalent; the closed forms were used here because
%   the algebra is short.
%
%   Degenerate interval: As in ssat.m, b == a leaves 1/(b-a) undefined, and
%   the convention
%       dssat(x,a,a) = 0
%   is adopted, since zeta(x;a,a) = sign(x-a) is piecewise constant and so
%   has derivative zero away from the single point a.
%
%   Inputs:
%       x      raw argument, any array shape
%       a, b   transition endpoints, scalars
%
%   Outputs:
%       dy_by_dx   d/dx of zeta(x;a,b), same shape as x
%
%   Called by DTHETA_STAR, DDTHETA_STAR, INNER_LOOP_CONTROLLER, and the local
%   JSSAT and JINTSSAT helpers of ALTERNATE_REFERENCE_FUNCTION.

if nargin  == 1
    s = x;
    ds_by_dx = 1;
elseif nargin == 2
    s = x./a;
    ds_by_dx = 1./a;
elseif nargin == 3
    s = (x-a)./(b-a);
    ds_by_dx = 1./(b-a);
    if a == b
        dy_by_dx = zeros(size(x));
        return
    end
else
    error('Bad number of arguments in dssat().')
end

dy_by_ds = (s>0).*(s<1).*(...
    exp(1./s + 1./(1 - s))./((exp(1./(1 - s)) + exp(1./s)).^2.*(1 - s).^2.*s.^2) + (2.*exp(1./s + 1./(1 - s)))./((exp(1./(1 - s)) + exp(1./s)).^2.*(1 - s).^2) - (2.*exp(1./s + 1./(1 - s)))./((exp(1./(1 - s)) + exp(1./s)).^2.*(1 - s).^2.*s)...
    );
dy_by_ds(s<=0) = 0;
dy_by_ds(s>=1) = 0;

dy_by_dx = dy_by_ds.*ds_by_dx;

dy_by_dx(x==-Inf) = 0;
dy_by_dx(x==Inf) = 0;

if any(isnan(dy_by_dx))
    dy_by_dx(isnan(dy_by_dx)) = 0;
end

end


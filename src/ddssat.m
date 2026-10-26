function ddy_by_ddx = ddssat(x,a,b)
%DDSSAT  Second derivative of ssat(x,a,b) with respect to x.
%   Companion to ssat.m, dssat.m, dddssat.m, ddddssat.m. Same calling
%   conventions:
%       ddssat(x)      ->  s = x,             ds/dx = 1
%       ddssat(x,a)    ->  s = x./a,          ds/dx = 1/a
%       ddssat(x,a,b)  ->  s = (x-a)./(b-a),  ds/dx = 1/(b-a)
%
%   IMPLEMENTATION NOTE
%   The inner step y(s) = phi./(phi+psi) is a reparameterized logistic,
%       y = 1/(1 + exp(w)),   w(s) = 1/s - 1/(1-s),
%   so the s-derivatives of y are logistic polynomials in (p,q) = (y,1-y)
%   times derivatives of w(s). To second order,
%       y'   = L1 w'
%       y''  = L2 w'^2 + L1 w''
%   with L1 = -pq and L2 = pq(1-2p). This file states dy/ds and ddy/dds as
%   single rational closed forms in the exponentials instead of forming the
%   polynomials explicitly; the two are equivalent, and the closed forms are
%   used here because the algebra is short. See DDDSSAT and DDDDSSAT for the
%   polynomial form.
%
%   Since s(x) is affine, all inner derivatives d^k s/dx^k vanish for k >= 2,
%   and the Faa di Bruno expansion of the second derivative retains only
%       ddy/dx^2 = (ddy/ds^2)(ds/dx)^2.
%   The p and q used in the polynomial form are built without subtraction
%   (q = psi/(phi+psi)) so they remain accurate in the tails where 1-p would
%   otherwise underflow. Both dy/ds and ddy/dds are zeroed outside the
%   support 0 < s < 1.
%
%   Degenerate interval: As in ssat.m and dssat.m, b == a leaves 1/(b-a)
%   undefined, and the convention
%       ddssat(x,a,a) = 0
%   is adopted, since the first derivative dssat(x;a,a) is zero away from the
%   single point a and so its derivative vanishes almost everywhere.
%
%   Inputs:
%       x      raw argument, any array shape
%       a, b   transition endpoints, scalars
%
%   Outputs:
%       ddy_by_ddx   d^2/dx^2 of zeta(x;a,b), same shape as x
%
%   Called by DDTHETA_STAR, and by the local JSSAT and JINTSSAT helpers of
%   ALTERNATE_REFERENCE_FUNCTION. Reached from MULTIROTOR_JUMP_GUARD and
%   MULTIROTOR_LYAPUNOV_FUNCTION indirectly, through DDTHETA_STAR.

if nargin  == 1
    s = x;
    ds_by_dx = 1;
    dds_by_ddx = 0;
elseif nargin == 2
    s = x./a;
    ds_by_dx = 1./a;
    dds_by_ddx = 0;
elseif nargin == 3
    s = (x-a)./(b-a);
    ds_by_dx = 1./(b-a);
    dds_by_ddx = 0;
    if a == b
        ddy_by_ddx = zeros(size(x));
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

ddy_by_dds = (s>0).*(s<1).*(...
    -(exp(1./(s - s.^2)).*(exp(1./(1 - s)).*(4.*s.^5 - 6.*s.^4 + 4.*s.^3 - 2.*s + 1) + exp(1./s).*(4.*s.^5 - 14.*s.^4 + 20.*s.^3 - 16.*s.^2 + 6.*s - 1)))./((exp(1./(1 - s)) + exp(1./s)).^3.*(s - 1).^4.*s.^4)...
    );
ddy_by_dds(s<=0) = 0;
ddy_by_dds(s>=1) = 0;

ddy_by_ddx = dy_by_ds.*dds_by_ddx + ddy_by_dds.*(ds_by_dx).^2;

ddy_by_ddx(x==-Inf) = 0;
ddy_by_ddx(x==Inf) = 0;

if any(isnan(ddy_by_ddx))
    ddy_by_ddx(isnan(ddy_by_ddx)) = 0;
end

end


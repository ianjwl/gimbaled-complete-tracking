function dddy_by_dddx = dddssat(x,a,b)
%DDDSSAT  Third derivative of ssat(x,a,b) with respect to x.
%   Companion to ssat.m, dssat.m, ddssat.m.  Same calling conventions:
%       dddssat(x)      ->  s = x,            ds/dx = 1
%       dddssat(x,a)    ->  s = x./a,         ds/dx = 1/a
%       dddssat(x,a,b)  ->  s = (x-a)./(b-a), ds/dx = 1/(b-a)
%
%   IMPLEMENTATION NOTE
%   The inner step y(s) = phi/(phi+psi), phi = exp(-1/s), psi = exp(-1/(1-s)),
%   is a reparameterized logistic:  y = 1/(1+exp(w)) with  w(s) = 1/s - 1/(1-s).
%   Hence the s-derivatives of y are logistic polynomials in (p,q) = (y,1-y)
%   times derivatives of w(s):
%       y'   = L1 w'
%       y''  = L2 w'^2 + L1 w''
%       y''' = L3 w'^3 + 3 L2 w' w'' + L1 w'''
%   with  L1 = -pq,  L2 = pq(1-2p),  L3 = -pq(1-6p+6p^2),  pq = p*q.
%   p and q are formed WITHOUT subtraction (q = psi/(phi+psi)) so they stay
%   accurate in the tails where 1-p would otherwise underflow.
%   The x-derivative follows from Faa di Bruno; s(x) is affine, so the inner
%   derivatives d^k s/dx^k vanish for k >= 2 (kept explicitly for clarity).

if nargin == 1
    s = x;              ds_by_dx = 1;        dds_by_ddx = 0;  ddds_by_dddx = 0;
elseif nargin == 2
    s = x./a;           ds_by_dx = 1./a;     dds_by_ddx = 0;  ddds_by_dddx = 0;
elseif nargin == 3
    s = (x-a)./(b-a);   ds_by_dx = 1./(b-a); dds_by_ddx = 0;  ddds_by_dddx = 0;
    if a == b
        dddy_by_dddx = zeros(size(x));
        return
    end
else
    error('Bad number of arguments in dddssat().')
end

in = (s > 0) & (s < 1);                 % support of the transition

% --- step value p = y and complement q = 1-y (subtraction-free) ----------
phi = exp(-1./s)    .* (s > 0);
psi = exp(-1./(1-s)).* ((1-s) > 0);
den = phi + psi;
p   = phi./den;
q   = psi./den;
pq  = p.*q;

% --- derivatives of w(s) = 1/s - 1/(1-s) ---------------------------------
wp   = -1./s.^2 - 1./(1-s).^2;
wpp  =  2./s.^3 - 2./(1-s).^3;
wppp = -6./s.^4 - 6./(1-s).^4;

% --- logistic-in-w derivatives (polynomials in p) ------------------------
L1 = -pq;
L2 =  pq.*(1 - 2*p);
L3 = -pq.*(1 - 6*p + 6*p.^2);

% --- y-derivatives w.r.t. s, zeroed outside the support ------------------
dy_by_ds     = L1.*wp;
ddy_by_dds   = L2.*wp.^2 + L1.*wpp;
dddy_by_ddds = L3.*wp.^3 + 3*L2.*wp.*wpp + L1.*wppp;

dy_by_ds(~in)     = 0;
ddy_by_dds(~in)   = 0;
dddy_by_ddds(~in) = 0;

% --- Faa di Bruno: third derivative w.r.t. x -----------------------------
dddy_by_dddx = dy_by_ds.*ddds_by_dddx ...
             + 3*ddy_by_dds.*ds_by_dx.*dds_by_ddx ...
             + dddy_by_ddds.*(ds_by_dx).^3;

dddy_by_dddx(x==-Inf) = 0;
dddy_by_dddx(x==Inf) = 0;

if any(isnan(dddy_by_dddx))
    % disp('dddy_by_dddx is NaN!')
    dddy_by_dddx(isnan(dddy_by_dddx)) = 0;
end

end
function ddddy_by_ddddx = ddddssat(x,a,b)
%DDDDSSAT  Fourth derivative of ssat(x,a,b) with respect to x.
%   Companion to ssat.m, dssat.m, ddssat.m, dddssat.m.  Same conventions:
%       ddddssat(x)      ->  s = x,            ds/dx = 1
%       ddddssat(x,a)    ->  s = x./a,         ds/dx = 1/a
%       ddddssat(x,a,b)  ->  s = (x-a)./(b-a), ds/dx = 1/(b-a)
%
%   IMPLEMENTATION NOTE
%   With y(s) = 1/(1+exp(w)), w(s) = 1/s - 1/(1-s), the s-derivatives of y are
%   logistic polynomials in (p,q) = (y,1-y) times derivatives of w(s):
%       y'    = L1 w'
%       y''   = L2 w'^2 + L1 w''
%       y'''  = L3 w'^3 + 3 L2 w' w'' + L1 w'''
%       y'''' = L4 w'^4 + 6 L3 w'^2 w'' + L2 (3 w''^2 + 4 w' w''') + L1 w''''
%   with  L1 = -pq,  L2 = pq(1-2p),  L3 = -pq(1-6p+6p^2),
%         L4 = pq(1-14p+36p^2-24p^3),  pq = p*q.
%   p and q are formed without subtraction (q = psi/(phi+psi)) for tail
%   accuracy.  s(x) is affine, so d^k s/dx^k = 0 for k >= 2 (kept explicitly).

if nargin == 1
    s = x;              ds_by_dx = 1;        dds_by_ddx = 0;  ddds_by_dddx = 0;  dddds_by_ddddx = 0;
elseif nargin == 2
    s = x./a;           ds_by_dx = 1./a;     dds_by_ddx = 0;  ddds_by_dddx = 0;  dddds_by_ddddx = 0;
elseif nargin == 3
    s = (x-a)./(b-a);   ds_by_dx = 1./(b-a); dds_by_ddx = 0;  ddds_by_dddx = 0;  dddds_by_ddddx = 0;
    if a == b
        ddddy_by_ddddx = zeros(size(x));
        return
    end
else
    error('Bad number of arguments in ddddssat().')
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
wp    =  -1./s.^2 -  1./(1-s).^2;
wpp   =   2./s.^3 -  2./(1-s).^3;
wppp  =  -6./s.^4 -  6./(1-s).^4;
wpppp =  24./s.^5 - 24./(1-s).^5;

% --- logistic-in-w derivatives (polynomials in p) ------------------------
L1 = -pq;
L2 =  pq.*(1 -  2*p);
L3 = -pq.*(1 -  6*p +  6*p.^2);
L4 =  pq.*(1 - 14*p + 36*p.^2 - 24*p.^3);

% --- y-derivatives w.r.t. s, zeroed outside the support ------------------
dy_by_ds       = L1.*wp;
ddy_by_dds     = L2.*wp.^2 + L1.*wpp;
dddy_by_ddds   = L3.*wp.^3 + 3*L2.*wp.*wpp + L1.*wppp;
ddddy_by_dddds = L4.*wp.^4 + 6*L3.*wp.^2.*wpp + L2.*(3*wpp.^2 + 4*wp.*wppp) + L1.*wpppp;

dy_by_ds(~in)       = 0;
ddy_by_dds(~in)     = 0;
dddy_by_ddds(~in)   = 0;
ddddy_by_dddds(~in) = 0;

% --- Faa di Bruno: fourth derivative w.r.t. x ----------------------------
ddddy_by_ddddx = dy_by_ds.*dddds_by_ddddx ...
               + ddy_by_dds.*(3*dds_by_ddx.^2 + 4*ds_by_dx.*ddds_by_dddx) ...
               + 6*dddy_by_ddds.*(ds_by_dx).^2.*dds_by_ddx ...
               + ddddy_by_dddds.*(ds_by_dx).^4;

ddddy_by_ddddx(x==-Inf) = 0;
ddddy_by_ddddx(x==Inf) = 0;

if any(isnan(ddddy_by_ddddx))
    % disp('ddddy_by_ddddx is NaN!')
    ddddy_by_ddddx(isnan(ddddy_by_ddddx)) = 0;
end

end
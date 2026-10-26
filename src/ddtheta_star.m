function [ddvalue,dvalue,value] = ddtheta_star(alpha,rho)
%DDTHETA_STAR  d^2/dalpha^2 of theta_star(alpha,rho), elementwise in alpha.
%   Also returns the first derivative and the value, which are computed
%   anyway. The inner loop needs ddt2_theta_star =
%   grad2_theta_star*ddt_alpha^2 + grad_theta_star*ddt2_alpha to form
%   ddt2_n_star, and thence the ddt_omega_star terms.
%
%   Construction. Same pieces as DTHETA_STAR, differentiated one order
%   further, so the value and the first derivative returned here agree
%   identically with THETA_STAR_FUNCTION and DTHETA_STAR:
%       I  = Z(alpha; alpha_pi1, alpha_plus)
%       D  = alpha - I, pinned to base for alpha >= alpha_plus
%       vartheta = 2*atan2(sqrt(D.*(1-D)), sqrt(I.*(1+I)))
%   The outer blend gives
%       value   = vartheta*sat + 2*pi*(1-sat)
%       dvalue  = vartheta'*sat + (vartheta - 2*pi)*dsat
%       ddvalue = vartheta''*sat + 2*vartheta'*dsat + (vartheta - 2*pi)*ddsat.
%
%   The two derivatives of vartheta are obtained from the cot(vartheta/2)
%   route used in DTHETA_STAR. With
%       g  = k*h,   k = sqrt(I),   h = sqrt((1+I)./Q),   Q = D.*(1-D),
%       L  = h'/h = 0.5*(v./(1+I) - dQ./Q),
%   the chain rule gives
%       g'  = h*(k' + k*L)
%       g'' = h*(k'' + 2*k'*L + k*(L^2 + L')),
%       vartheta'   = -2*g'/(1+g^2)
%       vartheta''  = -2*( g''/(1+g^2) - 2*g*(g')^2/(1+g^2)^2 ).
%   The quantities dk = v./(2*k), ddk = (dv/2 - dk^2)./k, dQ, and ddQ are
%   formed to avoid ever differencing two nearly equal numbers.
%
%   Domain. rho = U1/r is accepted on [0,1), as in DTHETA_STAR; rho = 0 is
%   the U_1 = 0 case and leaves all four breakpoints distinct.
%
%   Inputs:
%       alpha   normalized actuation norm ||a*||/r, any array shape
%       rho     U1/r, scalar in [0,1)
%
%   Outputs:
%       ddvalue  d^2/dalpha^2 of theta_star(alpha,rho), same shape as alpha
%       dvalue   d/dalpha of theta_star(alpha,rho)
%       value    theta_star(alpha,rho)
%
%   Requires: intssat, ssat, dssat, ddssat.
%
%   Called by INNER_LOOP_CONTROLLER.

if (rho >= 1)||(rho < 0); error('rho needs to be between 0 and 1!'); end
%   Also returns the first derivative and the value, which are computed anyway.
if (rho >= 1)||(rho < 0); error('rho needs to be between 0 and 1!'); end
base = (1-rho)/2;
alpha_2pi = base/4;
alpha_pi0 = base/2;
alpha_pi1 = 3/4*base;
alpha_plus = 5/4*base;

% inner saturation on [alpha_pi1, alpha_plus]
I  = intssat(alpha, alpha_pi1, alpha_plus);
v  = ssat(alpha, alpha_pi1, alpha_plus);
dv = dssat(alpha, alpha_pi1, alpha_plus);
D  = alpha - I;                     % s - lstarplus
D(alpha >= alpha_plus) = base;      % exact there

vartheta = 2*atan2(sqrt(D.*(1 - D)), sqrt(I.*(1 + I)));
vartheta(alpha <= alpha_pi1) = pi;  % degenerate triangle, including alpha = 0

% vartheta == pi identically where I == 0, so its derivatives are zero there
dvartheta  = zeros(size(alpha));
ddvartheta = zeros(size(alpha));
m = (alpha > alpha_pi1) & (I > 0);
if any(m(:))
    Im = I(m); vm = v(m); dvm = dv(m); Dm = D(m);
    k   = sqrt(Im);                     % sqrt(I)
    dk  = vm./(2*k);                    % ratio of two small, accurately computed numbers
    ddk = (dvm/2 - dk.^2)./k;
    Q   = Dm.*(1 - Dm);                 % >= alpha_pi1*(1 - alpha_pi1) > 0
    dQ  = (1 - vm).*(1 - 2*Dm);
    ddQ = -dvm.*(1 - 2*Dm) - 2*(1 - vm).^2;
    L   = 0.5*(vm./(1 + Im) - dQ./Q);   % h'/h
    dL  = 0.5*(dvm./(1 + Im) - (vm./(1 + Im)).^2 - ddQ./Q + (dQ./Q).^2);
    h   = sqrt((1 + Im)./Q);
    g   = k.*h;                         % cot(vartheta/2)
    dg  = h.*(dk + k.*L);
    ddg = h.*(ddk + 2*dk.*L + k.*(L.^2 + dL));
    G   = 1 + g.^2;
    dvartheta(m)  = -2*dg./G;
    ddvartheta(m) = -2*(ddg./G - 2*g.*dg.^2./G.^2);
end

% outer blend
sat   = ssat(alpha, alpha_2pi, alpha_pi0);
dsat  = dssat(alpha, alpha_2pi, alpha_pi0);
ddsat = ddssat(alpha, alpha_2pi, alpha_pi0);
value   = vartheta.*sat + 2*pi*(1 - sat);
dvalue  = dvartheta.*sat + (vartheta - 2*pi).*dsat;
ddvalue = ddvartheta.*sat + 2*dvartheta.*dsat + (vartheta - 2*pi).*ddsat;
end
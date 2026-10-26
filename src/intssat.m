function Z = intssat(x,a,b)
%INTSSAT  Integral of ssat(.,a,b) from a to x, elementwise in x (any array shape).
%   No approximation beyond adaptive quadrature at RelTol 1e-12. Exact
%   substitutions (tau = 1/w, w = 1/v + nu) and the reflection identity
%       int_0^u zeta0 = u - 1/2 + int_0^(1-u) zeta0,   u > 1/2,
%   reduce every element to one O(1), exponentially decaying integrand, so the
%   result has full relative accuracy everywhere, including just above x = a.
if a == b
    Z = max(x-b,0);
    return
end
u  = (x - a)./(b - a);
v  = min(u, 1 - u);                  % distance to the nearer end of the transition
ok = v >= 1/800;                     % below this the exact value underflows to 0
vc = ok.*min(v, 0.5) + (~ok)*0.5;    % benign placeholder for masked elements
Jt = zeros(size(vc));
for i = find(ok(:)).'               % one scalar adaptive integral per element in the transition
    Jt(i) = integral(@(nu) integrand(nu, vc(i)), 0, Inf, 'RelTol', 1e-12, 'AbsTol', 0);
end
F  = vc.^2 .* exp(-1./vc) .* Jt .* ok;   % int_0^v zeta0
Zn = (u <= 0.5).*F + (u > 0.5).*(u - 0.5 + F);     % also gives 0 for u <= 0 and u - 1/2 for u >= 1
Z  = (b - a).*Zn;
Z(x==-Inf) = 0;
Z(x==Inf) = Inf;
end

function y = integrand(nu, v)
% exp(-nu) * exp(w/(w-1)) / ((1 + exp(-X)) * (v*w)^2),  w = 1/v + nu,  X = w - w/(w-1)
w = 1./v + nu;
X = w - w./(w - 1);
y = exp(-nu + w./(w - 1)) ./ ((1 + exp(-X)) .* (v.*w).^2);
end

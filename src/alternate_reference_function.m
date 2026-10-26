function [x_user,dx_user,ddx_user,dddx_user,ddddx_user,not_free_falling] = alternate_reference_function(t)
%ALTERNATE_REFERENCE_FUNCTION  Reference position and its 1st-4th derivatives.
%   [x_user,dx_user,ddx_user,dddx_user,ddddx_user,not_free_falling] = ...
%       ALTERNATE_REFERENCE_FUNCTION(t)
%   returns the NED reference position x_user = [Y; X; -Z] and its exact time
%   derivatives up to fourth order, evaluated at time(s) t:
%
%       x_user      position            [Y;  X;  -Z]
%       dx_user     velocity            d/dt
%       ddx_user    acceleration        d2/dt2
%       dddx_user   jerk                d3/dt3
%       ddddx_user  snap                d4/dt4
%       not_free_falling   1 when actuated, 0 in the ballistic window
%
%   t may be scalar or a length-N vector. Each of x_user..ddddx_user is 3xN
%   (3x1 for scalar t); not_free_falling is 1xN.
%
%   Matches the revised Desmos "TAC2026 Flight Path". Derivatives are EXACT,
%   not finite-differenced: every factor is carried as a "jet"
%   [f; f'; f''; f'''; f''''] and combined with the Leibniz product rule
%   (see JMUL). Smooth ramps supply their derivatives via
%   ssat/dssat/ddssat/dddssat/ddddssat; the integrated ramp z_I via JINTSSAT
%   (d/dt z_I = ssat); polynomial factors in closed form.
%
%   REVISION NOTES (relative to the previous version):
%     * T5 = T4 + 2 (was +1), so T5=14, T6=17.5, T7=21.5, T8=23.5.
%     * New dwell offsets D_T3 = D_T6 = 0.5 introduce short holds:
%         - c_z descent uses [T3+D_T3, T4+D_T3]  -> brief hold at the climb top
%         - c_x settle and Y decay use the +D_T6 shifted edges
%     * X blend windows narrowed from +/-d to +/-d/2.
%
%   Requires on the MATLAB path:
%       ssat, dssat, ddssat, dddssat, ddddssat, intssat.

    t = t(:).';                 % row, so the 3-channel stacks are 3xN
    N = numel(t);

    %% --- Constants ------------------------------------------------------
    g   = 9.81;     % gravitational acceleration  [m/s^2]
    D   = 15;       % climb-amplitude parameter (Desmos D)
    d   = 1;        % width parameter of the smooth blend windows (Desmos d)
    DT3 = 0.5;      % dwell offset at the climb top   (Desmos D_{T3})
    DT6 = 0.5;      % dwell offset around T6          (Desmos D_{T6})

    %% --- Time knots (Desmos T0..T8) ------------------------------------
    T0 = 0;          % 0
    T1 = T0 + 2;     % 2
    T2 = T1 + 4;     % 6
    T3 = T2 + 3;     % 9
    T4 = T3 + 3;     % 12
    T5 = T4 + 2;     % 14
    T6 = T5 + 3.5;   % 17.5
    T7 = T6 + 4;     % 21.5
    T8 = T7 + 2;     % 23.5

    %% --- Lateral channel X(t) ------------------------------------------
    axJ = 3 * jssat(t, T1, T2);                     % 0 -> 3 ramp in
    cxJ = -3 * jssat(t, T7+DT6, T6+DT6);            % settle at -3 (shifted falling edge)
    kbx = -6 / (T6 - T5);                           % slope of bx
    bxJ        = zeros(5,N);                         % bx = 3 + kbx*(t-T5), affine
    bxJ(1,:)   = 3 + kbx*(t - T5);
    bxJ(2,:)   = kbx;

    XJ = jmul(axJ, jssat(t, T5+d/2, T5-d/2)) ...                                % t < T5
       + jmul(jmul(bxJ, jssat(t, T5-d/2, T5+d/2)), jssat(t, T6+d/2, T6-d/2)) ...% T5 < t < T6
       + jmul(cxJ, jssat(t, T6-d/2, T6+d/2));                                   % t > T6

    %% --- Lateral channel Y(t) ------------------------------------------
    YJ = 4 * jssat(t, T1, T2) - 4 * jssat(t, T6+DT6, T7+DT6);   % up to 4, back to 0

    %% --- Vertical channel Z(t):  Z = az + (bz-2) + (cz-(2+D)) + (dz-2) + (ez-2)
    azJ    =  2 * jssat(t, T0, T1);                 % initial rise to 2
    bz_m2J =  D * jssat(t, T2, T3);                 % (bz - 2):     climb  +D
    cz_mJ  = -D * jssat(t, T3+DT3, T4+DT3);         % (cz - (2+D)): descend -D (held D_T3)

    % ballistic arc:  dz - 2 = -(g/2) * z_I(.,T5-d/2,T5+d/2) * z_I(.,T6+d/2,T6-d/2)
    P = jintssat(t, T5 - d/2, T5 + d/2);            % rising  integrated ramp ~ max(0, t-T5)
    Q = jintssat(t, T6 + d/2, T6 - d/2);            % falling integrated ramp ~ min(0, t-T6)
    dz_m2J = -0.5 * g * jmul(P, Q);                 % (dz - 2)

    ez_m2J = -2 * jssat(t, T7, T8);                 % (ez - 2):     settle 2 -> 0

    ZJ = azJ + bz_m2J + cz_mJ + dz_m2J + ez_m2J;

    %% --- Assemble NED outputs:  x_user = [Y; X; -Z] at each order ------
    x_user     = [ YJ(1,:); XJ(1,:); -ZJ(1,:) ];
    dx_user    = [ YJ(2,:); XJ(2,:); -ZJ(2,:) ];
    ddx_user   = [ YJ(3,:); XJ(3,:); -ZJ(3,:) ];
    dddx_user  = [ YJ(4,:); XJ(4,:); -ZJ(4,:) ];
    ddddx_user = [ YJ(5,:); XJ(5,:); -ZJ(5,:) ];

    %% --- 0 if free-falling, 1 otherwise (value only) -------------------
    not_free_falling = ssat(t, T5+d/2, T5-d/2) + ssat(t, T6-d/2, T6+d/2);
end

% ========================================================================
function J = jssat(t,a,b)
%JSSAT  Jet of ssat(t,a,b): [value; d/dt; d2/dt2; d3/dt3; d4/dt4]  (5xN).
    J = [ ssat(t,a,b);
          dssat(t,a,b);
          ddssat(t,a,b);
          dddssat(t,a,b);
          ddddssat(t,a,b) ];
end

% ========================================================================
function J = jintssat(t,a,b)
%JINTSSAT  Jet of z_I(t,a,b) = int_a^t ssat(s,a,b) ds  (5xN).
%   Derivatives ARE the ssat ladder (fundamental theorem of calculus):
%       z_I' = ssat,  z_I'' = dssat,  z_I''' = ddssat,  z_I'''' = dddssat.
%   intssat.m requires its 2nd arg < 3rd arg.  When a > b (e.g.
%   z_I(t,T6+d/2,T6-d/2)) we use ssat(s,a,b) = 1 - ssat(s,b,a), giving
%       z_I(t,a,b) = (t - (a+b)/2) - intssat(t,b,a).
%   The derivative rows use the original (a,b), preserving the falling sign.
    if b > a
        I0 = intssat(t, a, b);
    else
        I0 = (t - 0.5*(a + b)) - intssat(t, b, a);
    end
    J = [ I0;
          ssat(t,a,b);
          dssat(t,a,b);
          ddssat(t,a,b);
          dddssat(t,a,b) ];
end

% ========================================================================
function w = jmul(u,v)
%JMUL  Leibniz product of two jets u,v (each 5xN) to 4th order.
%   (uv)^(n) = sum_{k=0}^{n} nchoosek(n,k) u^(k) v^(n-k)
    w = zeros(size(u));
    w(1,:) = u(1,:).*v(1,:);
    w(2,:) = u(2,:).*v(1,:) +   u(1,:).*v(2,:);
    w(3,:) = u(3,:).*v(1,:) + 2*u(2,:).*v(2,:) +   u(1,:).*v(3,:);
    w(4,:) = u(4,:).*v(1,:) + 3*u(3,:).*v(2,:) + 3*u(2,:).*v(3,:) +   u(1,:).*v(4,:);
    w(5,:) = u(5,:).*v(1,:) + 4*u(4,:).*v(2,:) + 6*u(3,:).*v(3,:) + 4*u(2,:).*v(4,:) + u(1,:).*v(5,:);
end

% function [x_user,dx_user,ddx_user,dddx_user,ddddx_user,not_free_falling] = alternate_reference_function(t)
% %ALTERNATE_REFERENCE_FUNCTION  Reference position and its 1st-4th derivatives.
% %   [x_user,dx_user,ddx_user,dddx_user,ddddx_user,not_free_falling] = ...
% %       ALTERNATE_REFERENCE_FUNCTION(t)
% %   returns the NED reference position x_user = [Y; X; -Z] and its exact time
% %   derivatives up to fourth order, evaluated at time(s) t:
% %
% %       x_user      position            [Y;  X;  -Z]
% %       dx_user     velocity            d/dt
% %       ddx_user    acceleration        d2/dt2
% %       dddx_user   jerk                d3/dt3
% %       ddddx_user  snap                d4/dt4
% %       not_free_falling   1 when actuated, 0 in the ballistic window
% %
% %   t may be scalar or a length-N vector. Each of x_user..ddddx_user is 3xN
% %   (3x1 for scalar t); not_free_falling is 1xN.
% %
% %   This matches the revised Desmos "TAC2026 Flight Path". Derivatives are
% %   EXACT, not finite-differenced: every factor is carried as a "jet"
% %   [f; f'; f''; f'''; f''''] and combined with the Leibniz product rule
% %   (see JMUL). Smooth-ramp factors supply their derivatives through
% %   ssat/dssat/ddssat/dddssat/ddddssat; polynomial factors supply theirs in
% %   closed form.
% %
% %   HEIGHT PROFILE (revised). The ballistic arc is now built from integrated
% %   ramps z_I(t,a,b) = int_a^t ssat(s,a,b) ds instead of an explicit
% %   parabola:
% %       dz - 2 = -(g/2) * z_I(t, T5-d/2, T5+d/2) * z_I(t, T6+d/2, T6-d/2)
% %   This still gives dz'' = -g (free fall) on the arc, but with C^inf
% %   entry/exit. Because d/dt z_I = ssat, the z_I jet reuses the ssat ladder
% %   for its derivatives; see JINTSSAT (which calls intssat for the value).
% %
% %   Requires on the MATLAB path:
% %       ssat, dssat, ddssat, dddssat, ddddssat, intssat.
% 
%     t = t(:).';                 % row, so the 3-channel stacks are 3xN
%     N = numel(t);
% 
%     %% --- Constants ------------------------------------------------------
%     g = 9.81;       % gravitational acceleration  [m/s^2]
%     D = 15;         % climb-amplitude parameter (Desmos D)
%     d = 1;          % half-width of the smooth blend windows (Desmos d)
% 
%     %% --- Time knots (Desmos T0..T8) ------------------------------------
%     T0 = 0;          % 0
%     T1 = T0 + 2;     % 2
%     T2 = T1 + 4;     % 6
%     T3 = T2 + 3;     % 9
%     T4 = T3 + 3;     % 12
%     T5 = T4 + 1;     % 13
%     T6 = T5 + 3.5;   % 16.5   (Desmos value; see note below)
%     T7 = T6 + 4;     % 20.5
%     T8 = T7 + 2;     % 22.5
% 
%     %% --- Lateral channel X(t) ------------------------------------------
%     axJ = 3 * jssat(t, T1, T2);                 % 0 -> 3 ramp in
%     cxJ = -3 * jssat(t, T7, T6);                % settle at -3 (falling edge)
%     kbx = -6 / (T6 - T5);                       % slope of bx
%     bxJ        = zeros(5,N);                     % bx = 3 + kbx*(t-T5), affine
%     bxJ(1,:)   = 3 + kbx*(t - T5);
%     bxJ(2,:)   = kbx;
% 
%     XJ = jmul(axJ, jssat(t, T5+d, T5-d)) ...                            % t < T5
%        + jmul(jmul(bxJ, jssat(t, T5-d, T5+d)), jssat(t, T6+d, T6-d)) ...% T5 < t < T6
%        + jmul(cxJ, jssat(t, T6-d, T6+d));                              % t > T6
% 
%     %% --- Lateral channel Y(t) ------------------------------------------
%     YJ = 4 * jssat(t, T1, T2) - 4 * jssat(t, T6, T7);   % up to 4, back to 0
% 
%     %% --- Vertical channel Z(t):  Z = az + (bz-2) + (cz-(2+D)) + (dz-2) + (ez-2)
%     azJ    =  2 * jssat(t, T0, T1);             % initial rise to 2
%     bz_m2J =  D * jssat(t, T2, T3);             % (bz - 2):     climb  +D
%     cz_mJ  = -D * jssat(t, T3, T4);             % (cz - (2+D)): descend -D
% 
%     % ballistic arc:  dz - 2 = -(g/2) * z_I(.,T5-d/2,T5+d/2) * z_I(.,T6+d/2,T6-d/2)
%     P = jintssat(t, T5 - d/2, T5 + d/2);        % rising  integrated ramp ~ max(0, t-T5)
%     Q = jintssat(t, T6 + d/2, T6 - d/2);        % falling integrated ramp ~ min(0, t-T6)
%     dz_m2J = -0.5 * g * jmul(P, Q);             % (dz - 2)
% 
%     ez_m2J = -2 * jssat(t, T7, T8);             % (ez - 2):     settle 2 -> 0
% 
%     ZJ = azJ + bz_m2J + cz_mJ + dz_m2J + ez_m2J;
% 
%     %% --- Assemble NED outputs:  x_user = [Y; X; -Z] at each order ------
%     x_user     = [ YJ(1,:); XJ(1,:); -ZJ(1,:) ];
%     dx_user    = [ YJ(2,:); XJ(2,:); -ZJ(2,:) ];
%     ddx_user   = [ YJ(3,:); XJ(3,:); -ZJ(3,:) ];
%     dddx_user  = [ YJ(4,:); XJ(4,:); -ZJ(4,:) ];
%     ddddx_user = [ YJ(5,:); XJ(5,:); -ZJ(5,:) ];
% 
%     %% --- 0 if free-falling, 1 otherwise (value only) ------------------
%     not_free_falling = ssat(t, T5+d, T5-d) + ssat(t, T6-d, T6+d);
% end
% 
% % ========================================================================
% function J = jssat(t,a,b)
% %JSSAT  Jet of ssat(t,a,b): [value; d/dt; d2/dt2; d3/dt3; d4/dt4]  (5xN).
%     J = [ ssat(t,a,b);
%           dssat(t,a,b);
%           ddssat(t,a,b);
%           dddssat(t,a,b);
%           ddddssat(t,a,b) ];
% end
% 
% % ========================================================================
% function J = jintssat(t,a,b)
% %JINTSSAT  Jet of z_I(t,a,b) = int_a^t ssat(s,a,b) ds  (5xN).
% %   By the fundamental theorem of calculus the derivatives ARE the ssat
% %   ladder:  z_I' = ssat,  z_I'' = dssat,  z_I''' = ddssat,  z_I'''' = dddssat.
% %
% %   intssat.m requires its 2nd arg < 3rd arg.  When a > b (a falling edge,
% %   e.g. z_I(t,T6+d/2,T6-d/2)) we use ssat(s,a,b) = 1 - ssat(s,b,a), giving
% %       z_I(t,a,b) = (t - (a+b)/2) - intssat(t,b,a).
% %   The derivative rows use the original (a,b) so the falling-edge sign is
% %   preserved automatically.
%     if b > a
%         I0 = intssat(t, a, b);
%     else
%         I0 = (t - 0.5*(a + b)) - intssat(t, b, a);
%     end
%     J = [ I0;
%           ssat(t,a,b);
%           dssat(t,a,b);
%           ddssat(t,a,b);
%           dddssat(t,a,b) ];
% end
% 
% % ========================================================================
% function w = jmul(u,v)
% %JMUL  Leibniz product of two jets u,v (each 5xN) to 4th order.
% %   (uv)^(n) = sum_{k=0}^{n} nchoosek(n,k) u^(k) v^(n-k)
%     w = zeros(size(u));
%     w(1,:) = u(1,:).*v(1,:);
%     w(2,:) = u(2,:).*v(1,:) +   u(1,:).*v(2,:);
%     w(3,:) = u(3,:).*v(1,:) + 2*u(2,:).*v(2,:) +   u(1,:).*v(3,:);
%     w(4,:) = u(4,:).*v(1,:) + 3*u(3,:).*v(2,:) + 3*u(2,:).*v(3,:) +   u(1,:).*v(4,:);
%     w(5,:) = u(5,:).*v(1,:) + 4*u(4,:).*v(2,:) + 6*u(3,:).*v(3,:) + 4*u(2,:).*v(4,:) + u(1,:).*v(5,:);
% end
% 
% % function [x_user,dx_user,ddx_user,dddx_user,ddddx_user,not_free_falling] = alternate_reference_function(t)
% % %POSITION_REFERENCE_DERIVATIVES  Reference position and its 1st-4th derivatives.
% % %   [x_user,dx_user,ddx_user,dddx_user,ddddx_user,not_free_falling] = ...
% % %       POSITION_REFERENCE_DERIVATIVES(t)
% % %   returns the NED reference position x_user = [Y; X; -Z] and its exact time
% % %   derivatives up to fourth order, evaluated at time(s) t:
% % %
% % %       x_user      position            [Y;  X;  -Z]
% % %       dx_user     velocity            d/dt
% % %       ddx_user    acceleration        d2/dt2
% % %       dddx_user   jerk                d3/dt3
% % %       ddddx_user  snap                d4/dt4
% % %       not_free_falling   1 when actuated, 0 in the ballistic window
% % %
% % %   t may be scalar or a length-N vector. Each of x_user..ddddx_user is 3xN
% % %   (3x1 for scalar t); not_free_falling is 1xN.
% % %
% % %   This mirrors position_reference.m exactly (same knots, gates, ENU->NED
% % %   swap). Derivatives are EXACT, not finite-differenced: every factor is
% % %   carried as a "jet" [f; f'; f''; f'''; f''''] and combined with the
% % %   Leibniz product rule (see JMUL below). The smooth-ramp factors supply
% % %   their own exact derivatives through ssat/dssat/ddssat/dddssat/ddddssat;
% % %   the polynomial factors (bx, dz) supply theirs in closed form.
% % %
% % %   Requires on the MATLAB path:
% % %       ssat, dssat, ddssat, dddssat, ddddssat.
% % 
% %     t = t(:).';                 % row, so the 3-channel stacks are 3xN
% %     N = numel(t);
% % 
% %     %% --- Constants ------------------------------------------------------
% %     g = 9.81;       % gravitational acceleration  [m/s^2]
% %     D = 15;         % climb-amplitude parameter (Desmos D)
% %     d = 1;          % half-width of the smooth blend windows (Desmos d)
% % 
% %     %% --- Time knots (Desmos T0..T8) ------------------------------------
% %     T0 = 0;         % 0
% %     T1 = T0 + 2;    % 2
% %     T2 = T1 + 4;    % 6
% %     T3 = T2 + 3;    % 9
% %     T4 = T3 + 3;    % 12
% %     T5 = T4 + 1;    % 13
% %     T6 = T5 + 6;    % 19
% %     T7 = T6 + 4;    % 22
% %     T8 = T7 + 2;    % 25
% % 
% %     %% --- Lateral channel X(t) ------------------------------------------
% %     axJ = 3 * jssat(t, T1, T2);                 % 0 -> 3 ramp in
% %     cxJ = -3 * jssat(t, T7, T6);                % settle at -3 (falling edge)
% %     kbx = -6 / (T6 - T5);                       % slope of bx
% %     bxJ        = zeros(5,N);                     % bx = 3 + kbx*(t-T5), affine
% %     bxJ(1,:)   = 3 + kbx*(t - T5);
% %     bxJ(2,:)   = kbx;
% % 
% %     XJ = jmul(axJ, jssat(t, T5+d, T5-d)) ...                            % t < T5
% %        + jmul(jmul(bxJ, jssat(t, T5-d, T5+d)), jssat(t, T6+d, T6-d)) ...% T5 < t < T6
% %        + jmul(cxJ, jssat(t, T6-d, T6+d));                              % t > T6
% % 
% %     %% --- Lateral channel Y(t) ------------------------------------------
% %     YJ = 4 * jssat(t, T1, T2) - 4 * jssat(t, T6, T7);   % up to 4, back to 0
% % 
% %     %% --- Vertical channel Z(t) (climb, then ballistic dip) -------------
% %     azJ    =  2 * jssat(t, T0, T1);             % initial rise to 2
% %     bz_m2J =  D * jssat(t, T2, T3);             % (bz - 2):       climb +D
% %     cz_mJ  = -D * jssat(t, T3, T4);             % (cz - (2+D)):   descend -D
% %     v5  = 0.5 * g * (T6 - T5);                  % launch rate of the ballistic arc
% %     dzJ        = zeros(5,N);                     % dz = 2 + v5*(t-T5) - g/2*(t-T5)^2
% %     dzJ(1,:)   = 2 + v5*(t - T5) - 0.5*g*(t - T5).^2;
% %     dzJ(2,:)   = v5 - g*(t - T5);
% %     dzJ(3,:)   = -g;
% %     ezJ    =  2 * jssat(t, T8, T7);             % settle back to 2 (falling edge)
% % 
% %     ZJ = azJ + bz_m2J + cz_mJ ...
% %        - 2 * jssat(t, T5, T5+2*d) ...                                       % remove baseline in free-fall
% %        + jmul(jmul(dzJ, jssat(t, T5, T5+2*d)), jssat(t, T6, T6-2*d)) ...    % ballistic arc, T5 < t < T6
% %        + jmul(ezJ, jssat(t, T6-2*d, T6));                                   % restore baseline after T6
% % 
% %     %% --- Assemble NED outputs:  x_user = [Y; X; -Z] at each order ------
% %     x_user     = [ YJ(1,:); XJ(1,:); -ZJ(1,:) ];
% %     dx_user    = [ YJ(2,:); XJ(2,:); -ZJ(2,:) ];
% %     ddx_user   = [ YJ(3,:); XJ(3,:); -ZJ(3,:) ];
% %     dddx_user  = [ YJ(4,:); XJ(4,:); -ZJ(4,:) ];
% %     ddddx_user = [ YJ(5,:); XJ(5,:); -ZJ(5,:) ];
% % 
% %     %% --- 0 if free-falling, 1 otherwise (value only) ------------------
% %     not_free_falling = ssat(t, T5+d, T5-d) + ssat(t, T6-d, T6+d);
% % end
% % 
% % % ========================================================================
% % function J = jssat(t,a,b)
% % %JSSAT  Jet of ssat(t,a,b): [value; d/dt; d2/dt2; d3/dt3; d4/dt4]  (5xN).
% %     J = [ ssat(t,a,b);
% %           dssat(t,a,b);
% %           ddssat(t,a,b);
% %           dddssat(t,a,b);
% %           ddddssat(t,a,b) ];
% % end
% % 
% % % ========================================================================
% % function w = jmul(u,v)
% % %JMUL  Leibniz product of two jets u,v (each 5xN) to 4th order.
% % %   (uv)^(n) = sum_{k=0}^{n} nchoosek(n,k) u^(k) v^(n-k)
% %     w = zeros(size(u));
% %     w(1,:) = u(1,:).*v(1,:);
% %     w(2,:) = u(2,:).*v(1,:) +   u(1,:).*v(2,:);
% %     w(3,:) = u(3,:).*v(1,:) + 2*u(2,:).*v(2,:) +   u(1,:).*v(3,:);
% %     w(4,:) = u(4,:).*v(1,:) + 3*u(3,:).*v(2,:) + 3*u(2,:).*v(3,:) +   u(1,:).*v(4,:);
% %     w(5,:) = u(5,:).*v(1,:) + 4*u(4,:).*v(2,:) + 6*u(3,:).*v(3,:) + 4*u(2,:).*v(4,:) + u(1,:).*v(5,:);
% % end
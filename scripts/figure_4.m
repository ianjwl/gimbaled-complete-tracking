function figure_4(save_figure)
%% ========================================================================
%  figure_4.m -- Geometry of the admissible set C(alpha) and target
%  a*(t)/r -- cases (a)-(d) (Panels export as figure_4_a ... figure_4_d.)
%
%  Everything is drawn in radius-normalized coordinates (lengths divided by
%  the target-set radius r):
%     alpha = a/r          normalized distance to the target-set centre
%     rho   = U1/r         normalized radius of the origin disk
%     D(alpha) = alpha*e1 + B^2   normalized target disk (unit radius)
%     rho*B^2                     normalized available-thrust disk
%
%  Four cases (a)-(d) for the location of a*(t)/r relative to D(alpha) and
%  rho*B^2. Each case is its own figure window (no subplots) so the panels
%  can be assembled as subfigures in LaTeX; the (a)-(d) sub-captions are
%  added in LaTeX, not here.
%
%  Every panel is exported on an identical fixed canvas: the plot box has
%  the same size and position in inches in all four PDFs, regardless of
%  where the annotations sit. Include each PDF at its natural width
%  (scale = 1) so the 8 pt fonts match the paper's other figures.
%
%  External helpers (must be on the MATLAB path):
%     theta_bar(alpha, rho)    admissible half-angle
%     theta_star_function(alpha, rho)   full angular extent (capped at pi when drawn)
%     tacStyle, tacAxes        shared figure styling
% ========================================================================

S = tacStyle;                 % shared style + LaTeX (CM) interpreters

% ---- Geometry (normalized by r) -----------------------------------------
rho = 0.6;         % normalized radius of the origin disk rho*B^2

% ---- Data shared by every panel -----------------------------------------
theta = linspace(-pi, pi, 2000);
[X, Y, C, map, nBands] = origin_disk_field(rho, S);

% D(alpha) label for panels (b)-(d): leader from the disk edge at -20 deg
% to a left-aligned label just right of the disk.
phiD = -20*pi/180;


%% ===== Panel (a):  alpha in [ 1 + rho , inf ) ==========================
%  Top band: a*(t)/r (left), theta*(alpha) (right, lower row),
%            theta_bar(alpha) (middle, upper row).
%  Bottom band: rho, alpha, D(alpha), C(alpha), left to right.
alpha = 2.5;

fig = figure('Theme','light'); hold on
draw_origin_disk(X, Y, C, map, nBands, rho, theta, S);
draw_calD(alpha, theta, S);
draw_calC(alpha, rho, S);
draw_axes();
draw_a_star_arrow(alpha, 0.5, S);
draw_theta_bar(alpha, rho, S, [0.16 0.98]);
draw_theta_star_function(alpha, rho, 0.3, '$\theta^{\star}(\alpha)$', S, [0.06 0.6], ["center" "bottom"]);

% --- dimension rulers unique to panel (a) -------------------------------
% ruler for alpha (label pulled left to make room for D(alpha))
xA = alpha/2 - 0.45;
plot([0 0],         -[0.025 0.125], 'k', 'LineWidth', S.lwMain);
plot([alpha alpha], -[0.025 0.125], 'k', 'LineWidth', S.lwMain);
plot(linspace(0,alpha,100), -0.075*ones(1,100), 'k', 'LineWidth', S.lwMain);
text(xA, -0.7, '$\alpha$', 'Interpreter','latex', 'FontSize',S.fsAnnot, ...
     'HorizontalAlignment','center', 'VerticalAlignment','top');
plot([alpha/2 xA], [-0.075 -0.75], 'k', 'LineWidth', S.lwLead);
% ruler for rho
plot([0 0],     -[0.025 0.125], 'k', 'LineWidth', S.lwMain);
plot([-rho -rho], -[0.025 0.125], 'k', 'LineWidth', S.lwMain);
plot(linspace(-rho,0,50), -0.075*ones(1,50), 'k', 'LineWidth', S.lwMain);
text(-rho/2-0.15, -0.7, '$\rho$', 'Interpreter','latex', 'FontSize',S.fsAnnot, ...
     'HorizontalAlignment','center', 'VerticalAlignment','top');
plot([-rho/2 -rho/2-0.15], [-0.075 -0.7], 'k', 'LineWidth', S.lwLead);

draw_a_star_label(alpha, [0.75 1.1], "bottom", S);
label_calD(alpha, -125*pi/180, 0.9, [1.65 -1.05], "center", "top", S);
label_calC(alpha, rho, [2.75 -1.25], S);
finalize_axes(fig, S);
if save_figure
    exportFig(fig, fullfile('figures','figure_4_a'));
end

%% ===== Panel (b):  alpha in [ sqrt(1 + rho^2) , 1 + rho ) ==============
alpha = 1.45;

fig = figure('Theme','light'); hold on
draw_origin_disk(X, Y, C, map, nBands, rho, theta, S);
draw_calD(alpha, theta, S);
draw_calC(alpha, rho, S);
draw_axes();
draw_a_star_arrow(alpha, 2, S);
draw_theta_bar(alpha, rho, S);
draw_theta_star_function(alpha, rho, 0.25, '$\theta^{\star}(\alpha)$', S);

draw_a_star_label(alpha, [alpha/2-0.25, -0.925], "top", S);
label_calD(alpha, phiD, 0.9, [alpha+1.35 -0.4], "left", "middle", S);
label_calC(alpha, rho, [2.6 -1.05], S);
finalize_axes(fig, S);
if save_figure
    exportFig(fig, fullfile('figures','figure_4_b'));
end

%% ===== Panel (c):  alpha in ( 1 - rho , sqrt(1 + rho^2) ) =============
alpha = 0.6;

fig = figure('Theme','light'); hold on
draw_origin_disk(X, Y, C, map, nBands, rho, theta, S);
draw_calD(alpha, theta, S);
draw_calC(alpha, rho, S);
draw_axes();
draw_a_star_arrow(alpha, 2, S);
draw_theta_bar(alpha, rho, S);
draw_theta_star_function(alpha, rho, 0.2, '$\theta^{\star}(\alpha)$', S);

draw_a_star_label(alpha, [alpha/2-0.35, -0.925], "top", S);
label_calD(alpha, phiD, 0.85, [alpha+1.35 -0.4], "left", "middle", S);
label_calC(alpha, rho, [1.8 -1.05], S);
finalize_axes(fig, S);
if save_figure
    exportFig(fig, fullfile('figures','figure_4_c'));
end

%% ===== Panel (d):  alpha in [ 0 , 1 - rho ] ===========================
%  Every direction is admissible (theta_bar = Inf), so no theta_bar ray is
%  drawn. At alpha = 0.05 = alpha_2pi, theta_star = 2*pi > pi, so the min
%  is active and C(alpha) is a full annulus.
alpha = 0.05;

fig = figure('Theme','light'); hold on
draw_origin_disk(X, Y, C, map, nBands, rho, theta, S);
draw_calD(alpha, theta, S);
draw_calC(alpha, rho, S);
draw_axes();
draw_a_star_arrow(alpha, 2, S);
draw_theta_star_function(alpha, rho, 0.1, '${\rm min}\{\pi,\theta^{\star}(\alpha)\}$', S);

draw_a_star_label(alpha, [alpha/2-0.35, -1.0], "top", S);   % left, clear of the C(alpha) leader
label_calD(alpha, phiD, 0.25, [alpha+1.35 -0.4], "left", "middle", S);
label_calC(alpha, rho, [1.6 -1.05], S);
finalize_axes(fig, S);
if save_figure
    exportFig(fig, fullfile('figures','figure_4_d'));
end

end

%% ======================= Local drawing helpers ==========================

function exportFig(fig, name)
%EXPORTFIG  Vector PDF (for the paper) + high-res PNG (for quick preview).
%  'Padding','figure' exports the whole fixed canvas set in finalize_axes
%  instead of cropping to the content, so all panels share one bounding box.
    try
        exportgraphics(fig, [name '.pdf'], 'ContentType','vector', 'Padding','figure');
    catch
        % Releases without the Padding option: print the full canvas instead.
        warning('exportFig:noPadding', ...
            'exportgraphics has no Padding option here; falling back to print.');
        fig.PaperUnits    = 'inches';
        fig.PaperSize     = fig.Position(3:4);
        fig.PaperPosition = [0 0 fig.Position(3:4)];
        print(fig, name, '-dpdf', '-vector');
    end
end

function [X, Y, C, map, nBands] = origin_disk_field(rho, S)
%ORIGIN_DISK_FIELD  Striped-grey field used to shade the rho-disk.
    nPairs = 5;
    nBands = 2*nPairs;
    map    = repmat([S.grey1; S.grey2], nPairs, 1);

    ehat = [ cos(pi/4)  sin(pi/4)];
    nhat = [-sin(pi/4)  cos(pi/4)];
    s    = rho * linspace(-1, 1, 101)';

    gx = s*ehat;
    gy = s*nhat;
    X  = gx(:,1) + gy(:,1)';
    Y  = gx(:,2) + gy(:,2)';

    C = nhat(1)*X + nhat(2)*Y;
    C(X.^2 + Y.^2 > rho^2) = NaN;
    C = C / max(C(:));
end

function draw_origin_disk(X, Y, C, map, nBands, rho, theta, S)
%DRAW_ORIGIN_DISK  Striped-grey fill of rho*B^2 plus its black boundary.
    [~, h] = contourf(X, Y, C, -1:(2/nBands):1);
    colormap(map);
    h.EdgeAlpha = 0;
    plot(rho*cos(theta), rho*sin(theta), 'k', 'LineWidth', S.lwThin);
end

function draw_calD(alpha, theta, S)
%DRAW_CALD  Light-green fill of the unit target disk D(alpha), with boundary.
    xs = cos(theta) + alpha;
    ys = sin(theta);
    fill(xs, ys, S.green, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    plot(xs, ys, 'k', 'LineWidth', S.lwThin);
end

function [vth, rin, rout] = C_geometry(alpha, rho)
%C_GEOMETRY  Half-angle and inner/outer radii of the annular section C(alpha).
    vth  = min(pi, theta_star_function(alpha, rho));
    disc = max(0, 1 - alpha^2*sin(vth)^2);
    u1p  = alpha*cos(vth) + sqrt(disc);
    u1m  = alpha*cos(vth) - sqrt(disc);
    rin  = max(rho, u1m);
    rout = u1p;
end

function draw_calC(alpha, rho, S)
%DRAW_CALC  Solid-blue fill of the admissible set C(alpha) (annular section).
    [vth, rin, rout] = C_geometry(alpha, rho);
    phi = linspace(-vth, vth, 600);
    xv  = [rout*cos(phi), fliplr(rin*cos(phi))];
    yv  = [rout*sin(phi), fliplr(rin*sin(phi))];
    fill(xv, yv, S.blue, 'EdgeColor', 'k', 'FaceAlpha', 0.8);
end

function draw_axes()
%DRAW_AXES  Faint Cartesian axes through the origin.
    xl = xline(0, 'k'); xl.Alpha = 0.35;
    yl = yline(0, 'k'); yl.Alpha = 0.35;
end

function draw_a_star_arrow(alpha, headSize, S)
%DRAW_A_STAR_ARROW  Arrow 0 -> (alpha,0) for a*(t)/r; a dot at the origin if alpha = 0.
    if alpha > 0
        quiver(0, 0, alpha, 0, 'k', 'LineWidth', S.lwMain, ...
               'AutoScale', 'off', 'MaxHeadSize', headSize);
    else
        plot(0, 0, '.k', 'MarkerSize', 12);
    end
end

function draw_theta_bar(alpha, rho, S, off)
%DRAW_THETA_BAR  Admissible half-angle theta_bar(alpha): two rays + labelled arc.
%  OFF = [dx dy] places the label's bottom-centre relative to the arc
%  midpoint (default [0.1 0.7]).
    if nargin < 4, off = [0.1 0.7]; end
    tb = theta_bar(alpha, rho);
    if ~isfinite(tb), return; end

    rt  = alpha + [0.1 0.2];
    tkx = [1; 0]*rt;
    tkt = [cos(tb); sin(tb)]*rt;
    arc = (alpha + 0.15)*[cos(linspace(0,tb,500)'), sin(linspace(0,tb,500)')];

    plot([0 alpha*cos(tb)], [0 alpha*sin(tb)], 'k', 'LineWidth', S.lwMain);
    plot(tkx(1,:), tkx(2,:), 'k', 'LineWidth', S.lwMain);
    plot(tkt(1,:), tkt(2,:), 'k', 'LineWidth', S.lwMain);
    plot(arc(:,1), arc(:,2), 'k', 'LineWidth', S.lwMain);

    p0 = [cos(tb/2); sin(tb/2)]*(alpha + 0.15);   % arc midpoint
    lp = p0 + off(:);
    text(lp(1), lp(2), '$\bar{\theta}(\alpha)$', 'Interpreter','latex', ...
         'FontSize',S.fsAnnot, 'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom');
    plot([p0(1) lp(1)], [p0(2) lp(2)], 'k', 'LineWidth', S.lwLead);
end

function draw_theta_star_function(alpha, rho, c, lbl, S, off, align)
%DRAW_THETA_STAR  Mark min(pi, theta_star_function(alpha)) at C(alpha)'s outer edge,
%  labelled with the LaTeX string LBL. OFF = [dx dy] places the label
%  anchor relative to the arc midpoint (default [0.5 0.7]). With no ALIGN
%  the anchor is the label's top-left corner and the leader meets its left
%  edge; ALIGN = [halign valign] sets the alignment explicitly, and the
%  leader then stops just short of that anchor.
    if nargin < 6 || isempty(off), off = [0.5 0.7]; end
    [vth, ~, baseR] = C_geometry(alpha, rho);

    rt  = baseR + [c-0.05 c+0.05];
    tkx = [1; 0]*rt;
    tkt = [cos(vth); sin(vth)]*rt;
    arc = (baseR + c)*[cos(linspace(0,vth,500)'), sin(linspace(0,vth,500)')];

    plot([0 baseR*cos(vth)], [0 baseR*sin(vth)], 'k', 'LineWidth', S.lwMain);
    plot(tkx(1,:), tkx(2,:), 'k', 'LineWidth', S.lwMain);
    plot(tkt(1,:), tkt(2,:), 'k', 'LineWidth', S.lwMain);
    plot(arc(:,1), arc(:,2), 'k', 'LineWidth', S.lwMain);

    p0 = [cos(vth/2); sin(vth/2)]*(baseR + c);   % arc midpoint
    lp = p0 + off(:);
    if nargin < 7
        halign = "left";  valign = "top";
        le = [lp(1), lp(2)-0.3];
    else
        halign = align(1);  valign = align(2);
        le = leader_end(lp, halign, valign);
    end
    text(lp(1), lp(2), lbl, 'Interpreter','latex', ...
         'FontSize',S.fsAnnot, 'HorizontalAlignment',halign, ...
         'VerticalAlignment',valign);
    plot([p0(1) le(1)], [p0(2) le(2)], 'k', 'LineWidth', S.lwLead);
end

function draw_a_star_label(alpha, pos, valign, S)
%DRAW_A_STAR_LABEL  a*(t)/r text at POS with a thin leader to the arrow midpoint.
    text(pos(1), pos(2), '${\rm a}^{\star}(t)/r$', 'Interpreter','latex', ...
         'FontSize',S.fsAnnot, 'HorizontalAlignment','center', ...
         'VerticalAlignment',valign);
    plot([alpha/2 pos(1)], [0 pos(2)-0.1], 'k', 'LineWidth', S.lwLead);
end

function label_calD(alpha, phi, k, pos, halign, valign, S)
%LABEL_CALD  D(alpha) label at POS, with a leader from the disk boundary
%  point at angle PHI (measured from the disk centre).
    p0 = [alpha + k*cos(phi), k*sin(phi)];
    le = leader_end(pos, halign, valign);
    plot([p0(1) le(1)], [p0(2) le(2)], 'k', 'LineWidth', S.lwLead);
    text(pos(1), pos(2), '$\mathcal{D}(\alpha)$', 'Interpreter','latex', ...
         'HorizontalAlignment',halign, 'VerticalAlignment',valign, ...
         'FontSize',S.fsAnnot);
end

function label_calC(alpha, rho, pos1, S)
%LABEL_CALC  C(alpha) label with a leader into the annulus section.
    [vth, rin, rout] = C_geometry(alpha, rho);
    rmid = (rin+rout)/2;
    pos0 = rmid*[cos(vth/2), -sin(vth/2)];
    le   = leader_end(pos1, "left", "middle");
    text(pos1(1), pos1(2), '$\mathcal{C}(\alpha)$', 'Interpreter','latex', ...
         'FontSize',S.fsAnnot, 'HorizontalAlignment','left', ...
         'VerticalAlignment','middle');
    plot([pos0(1) le(1)], [pos0(2) le(2)], 'k', 'LineWidth', S.lwLead);
end

function le = leader_end(pos, halign, valign)
%LEADER_END  Point a small gap outside a label's anchor, so a leader line
%  stops short of the text instead of touching it.
    g  = 0.08;                     % gap in data units
    le = pos(:).';
    switch halign
        case "left",  le(1) = le(1) - g;
        case "right", le(1) = le(1) + g;
    end
    switch valign
        case "top",    le(2) = le(2) + g;
        case "bottom", le(2) = le(2) - g;
    end
end

function finalize_axes(fig, S)
%FINALIZE_AXES  Common limits, normalized ticks/labels, shared axis style,
%  and a fixed canvas: the plot box is pinned to the same inch position in
%  every panel, so the exported PDFs are identical in size and frame.
    ax = fig.CurrentAxes;
    axis(ax, 'equal');
    xlim(ax, [-1.5 4]);
    ylim(ax, [-1.5 2]);

    xtks = -1:4;  xticks(ax, xtks);
    ytks = -1:2;  yticks(ax, ytks);
    % math mode so negative ticks get a true minus sign, not a hyphen
    xticklabels(ax, arrayfun(@(x) sprintf('$%g$', x), xtks, 'UniformOutput', false));
    yticklabels(ax, arrayfun(@(x) sprintf('$%g$', x), ytks, 'UniformOutput', false));

    tacAxes(ax, S);

    % Margins [in] around the plot box: room for the y tick labels (left),
    % x tick labels (bottom), and the half-labels at the corners (top/right).
    mL = 0.28;  mB = 0.20;  mR = 0.05;  mT = 0.06;
    w  = S.panelW - mL - mR;
    h  = w*(3.5/5.5);              % matches the [-1.5,4] x [-1.5,2] data aspect
    fig.Units    = 'inches';
    fig.Position = [1 1 S.panelW h+mB+mT];
    ax.Units         = 'inches';
    ax.InnerPosition = [mL mB w h];
end


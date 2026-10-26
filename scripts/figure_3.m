function figure_3(save_figure)
%% ======================================================================
%  figure_3.m  --  admissible theta_bar(s) and bound theta*(s)
%  Restyled to match the TAC sim figures (single column, Computer Modern).
%
%  The alpha breakpoints (alpha_2pi, alpha_pi0, alpha_pi1, alpha_plus) and
%  sqrt(1+rho^2) keep their tick marks, but their labels are drawn as
%  on-plot text with leader lines to the corresponding points on theta*
%  (for the alphas) and theta_bar (for sqrt(1+rho^2)).
%
%  External helpers required on the path:
%     theta_bar(s, rho), theta_star_function(s, rho), tacStyle, tacAxes
%% ======================================================================
S = tacStyle;                       % shared style + LaTeX (CM) interpreters
rho = 0.6;
s01 = 1 - rho;
s12 = sqrt(1 + rho^2);

% --- theta_bar branch (with vertical asymptote drawn at s = s01) --------
x    = [s01, linspace(s01, s12, 1e3), linspace(s12, 2.1, 1e3)];
y    = theta_bar(x, rho);
y(1) = 1e3;                         % spike -> near-vertical asymptote
x1 = x;  y1 = y;

% --- theta_star branch --------------------------------------------------
x = linspace(0, 2.1, 1e4);
[y,alpha_2pi,alpha_pi0,alpha_pi1,alpha_plus] = theta_star_function(x, rho);
x2 = x;  y2 = y;

%% --- figure ------------------------------------------------------------
fig = figure('Theme','light','Units','inches','Position',[1 1 S.colW 1]);
ax  = axes(fig); hold(ax,'on');
plot(ax, x1, y1, 'Color', S.blue, 'LineWidth', S.lwMain, 'DisplayName','$\bar{\theta}(\alpha)$');
plot(ax, x2, y2, 'k',             'LineWidth', S.lwMain, 'DisplayName','$\theta^{\star}(\alpha)$');
hold(ax,'off');
xlim(ax,[0 2.1]);
ylim(ax,[0 5*pi/2]);

% --- x ticks: marks at every breakpoint, axis labels only for the plain ones
xtPlain = [0 1 2 1-rho];
xlPlain = ["$0$" "$1$" "$2$" "$1-\rho$"];

% annotated breakpoints: tick mark on the axis, label drawn as text()
xtA = [alpha_2pi alpha_pi0 alpha_pi1 alpha_plus s12];
xlA = ["$\alpha_{2\pi}$" "$\alpha_{\pi,0}$" "$\alpha_{\pi,1}$" "$\alpha_{+}$" ...
       "$\sqrt{1+\rho^{2}}$"];

[xt, ord] = sort([xtPlain xtA]);
xl = [xlPlain repmat("", size(xtA))];
xticks(ax, xt);  xticklabels(ax, xl(ord));  xtickangle(ax, 0);
yticks(ax, [0 pi 2*pi]);  yticklabels(ax, ["$0$" "$\pi$" "$2\pi$"]);
tacAxes(ax, S);

% --- alpha breakpoint labels as on-plot text with leader lines ----------
% Leaders run from each label to its breakpoint on the relevant curve:
% the alphas to (alpha_i, theta*(alpha_i)), and (1+rho^2)^(1/2) to the
% branch switch of theta_bar at (s12, theta_bar(s12)).
yA = [theta_star_function(xtA(1:4), rho)  theta_bar(s12, rho)];

% Label anchors (data units) -- tune by eye. The leader ends at the anchor;
% hAlign/vAlign say which edge/corner of the text sits on the anchor.
labPos = [0.18 7.05              % alpha_2pi  : top strip, above the pi labels
          0.17 4.2               % alpha_pi0  : above the pi plateau
          0.31 4.2               % alpha_pi1  : above the pi plateau
          0.17 1.7               % alpha_plus : below the falling branch
          1.05 2.3];             % s12        : above the theta_bar curve
hAlign = ["left"   "center" "center" "center" "center"];
vAlign = ["middle" "bottom" "bottom" "top"    "bottom"];

hold(ax,'on');
for i = 1:numel(xtA)
    plot(ax, [xtA(i) labPos(i,1)], [yA(i) labPos(i,2)], 'b', ...
        'LineWidth', S.lwLead/2, 'HandleVisibility', 'off');   % keep out of legend
    text(ax, labPos(i,1), labPos(i,2), xlA(i), 'Interpreter', 'latex', ...
        'FontSize', S.fsAnnot, ...
        'HorizontalAlignment', hAlign(i), 'VerticalAlignment', vAlign(i));
end
hold(ax,'off');

legend(ax, 'Location','northeast', 'FontSize', S.fsLeg);

%% --- export ------------------------------------------------------------
if save_figure
    exportgraphics(fig, fullfile('figures','figure_3.pdf'), 'ContentType','vector');
end
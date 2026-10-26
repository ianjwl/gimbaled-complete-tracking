function figure_2(save_figure)
%% ======================================================================
%  figure_2.m  --  setup schematic for the a*(t) geometry
%  Shows the commanded vector a(t), its angle theta(t), the resulting
%  command a*(t), and the a / U1 dimension rulers, relative to the target
%  set A = a + r*B^2 and the origin thrust disk U1*B^2.
%
%  Four panels (a)-(d) over a_vals; panel (a) carries the extra a(t),
%  theta(t), and ruler annotations. Each panel is its own figure window
%  (no subplots) for assembly as subfigures in LaTeX.
%
%  Styling centralized in tacStyle/tacAxes to match the paper's figures.
%  Each panel is sized for one cell of a 2x2 grid in a single column
%  (~half column wide, S.panelW). Include each PDF at that natural width
%  (scale = 1) so the 8 pt fonts match the paper's other figures.
%
%  External helpers on the path: tacStyle, tacAxes.
%% ======================================================================

S = tacStyle;                   % shared style + LaTeX (CM) interpreters

% ---- Parameters ---------------------------------------------------------
rho = 0.6;                      % origin circle radius
r = 1;                          % radius of calA
critical_alpha_vals = [1 - rho  sqrt(1 + rho^2)  1 + rho];
alpha_vals = [2.5 1.45 0.6 0.2];
u1byr = 1.8;                    % length of the illustrative a(t) (panel a)

theta = linspace(-pi, pi, 2000);

% ---- D(alpha) label placement, one row per panel -----------------------
D_phi = [-70 -60 -50 -60]*pi/180;   % leader attach angle on the disk boundary
D_pos = [3.1  -1.75                 % label anchor (top-centre of the text)
         2.5  -2.1
         1.9  -2.1
         1.3  -1.75];

% ---- Striped-grey field for the origin disk U1*B^2 ----------------------
pair_count = 5;
bar_count  = 2*pair_count;
map = repmat([S.grey1; S.grey2], [pair_count 1]);
dirlines = pi/4;
xlines = rho*linspace(-1,1,1e2+1)';
ylines = rho*linspace(-1,1,1e2+1)';
dirx = [cos(dirlines) sin(dirlines)];
diry = [-sin(dirlines) cos(dirlines)];
xlines_vecs = xlines*dirx;
ylines_vecs = ylines*diry;
xlines_vecs_x = xlines_vecs(:,1);
xlines_vecs_y = xlines_vecs(:,2);
ylines_vecs_x = ylines_vecs(:,1)';
ylines_vecs_y = ylines_vecs(:,2)';
X = xlines_vecs_x + ylines_vecs_x;
Y = xlines_vecs_y + ylines_vecs_y;
C = diry(1)*X + diry(2)*Y;
C(X.^2+Y.^2>rho^2) = nan;
C = C/max(C(:));

letter = 'abcd';

for k = 1:numel(alpha_vals)
    fig = figure('Theme','light');

    alpha = alpha_vals(k)*r;

    %% --- Origin disk U1*B^2 (striped grey) ---
    [~,cf] = contourf(X,Y,C,-1:(2/bar_count):1);
    colormap(map)
    cf.EdgeAlpha = 0;
    axis equal
    hold on

    ndash = 20;
    sstep = ceil(numel(theta)/ndash);
    iinds = 1:numel(theta);
    miinds = mod(iinds-1,2*sstep)+1;
    indss = miinds > sstep;
    zros = zeros(size(theta));
    zros(indss) = nan;
    plot(rho*cos(theta) + zros, rho*sin(theta) + zros, 'k', 'LineWidth', S.lwThin);

    %% --- Target disk A = a + r*B^2 ---
    xs = cos(theta) + alpha;
    ys = sin(theta);
    fill(xs, ys, S.green, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    plot(xs, ys, 'k', 'LineWidth', S.lwThin);

    %% --- Faint Cartesian axes ---
    xl = xline(0,'k'); xl.Alpha = 0.35;
    yl = yline(0,'k'); yl.Alpha = 0.35;

    %% --- a*(t) arrow (or dot at the origin) ---
    if alpha > 0
        if k == 1; mk = 0.5; else; mk = 2; end
        quiver(0,0,alpha,0,'k','LineWidth',S.lwMain,'AutoScale','off','MaxHeadSize',mk)
    else
        plot(0,0,'.k','MarkerSize',8)
    end

    %% --- theta_bar rays + labelled arc (panels a-c) ---
    if k <= 3
        if alpha >= sqrt(1 + rho^2)
            theta_bar = asin(1/alpha);
        elseif alpha > 1 - rho
            theta_bar = acos((alpha^2 + rho^2 - 1)/(2*rho*alpha));
        else
            theta_bar = nan; % all angles admissible
        end

        b1 = [1; 0]*(alpha + [0.1 0.2]);
        b2 = [cos(theta_bar); sin(theta_bar)]*(alpha + [0.1 0.2]);
        b3 = (alpha + 0.15)*[cos(linspace(0,theta_bar,500)') sin(linspace(0,theta_bar,500)')];
        plot([0 alpha*cos(theta_bar)],[0 alpha*sin(theta_bar)], 'k', 'LineWidth', S.lwMain)
        plot(b1(1,:), b1(2,:), 'k', 'LineWidth', S.lwMain)
        plot(b2(1,:), b2(2,:), 'k', 'LineWidth', S.lwMain)
        plot(b3(:,1), b3(:,2), 'k', 'LineWidth', S.lwMain)

        label_point = [cos(theta_bar/2); sin(theta_bar/2)]*(alpha + 0.15) + [0.1;1];
        text(label_point(1), label_point(2), ...
            '$\bar{\theta}(\alpha)$', ...
            'Interpreter','latex', 'FontSize',S.fsAnnot, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle');
        plot([label_point(1)-0.1 label_point(1)],[label_point(2)-1 label_point(2)-0.4],'k','LineWidth',S.lwLead)
    end

    %% --- Extra annotations unique to panel (a) ---
    if k == 1
        % commanded vector a(t)
        q = 2*pi/3;
        quiver(0,0,u1byr*cos(q),u1byr*sin(q),'k','LineWidth',S.lwMain,'AutoScale','off','MaxHeadSize',0.5)
        text(1.1*u1byr*cos(q)+0.15,1.3*u1byr*sin(q)+0.1, ...
                '${\rm a}(t)/r$', ...
                'Interpreter','latex', 'FontSize',S.fsAnnot, ...
                'HorizontalAlignment','right', 'VerticalAlignment','middle');
        plot(u1byr*cos(q),u1byr*sin(q),'.k','MarkerSize',8)

        % angle theta(t)
        cq = [0.1 0.1];
        c1 = cq + [1; 0]*(0.7 + [-0.1 0.1]);
        c2 = cq + [cos(0.9*q); sin(0.9*q)]*(0.7 + [-0.1 0.1]);
        c2_int = cq + [cos(q); sin(q)]*(0.7);
        c2 = [c2(:,1)  c2_int c2(:,2)];
        c3 = cq + (0.7)*[cos(linspace(0,q,400)') sin(linspace(0,q,400)')];
        plot(c1(1,:), c1(2,:), 'k', 'LineWidth', S.lwMain)
        plot(c2(1,:), c2(2,:), 'k', 'LineWidth', S.lwMain)
        plot(c3(:,1), c3(:,2), 'k', 'LineWidth', S.lwMain)
        text(cos(q/2)*1.3+0.1, sin(q/2)*1.3, ...
            '$\theta(t)$', ...
            'Interpreter','latex', 'FontSize',S.fsAnnot, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle');

        % dimension ruler for a
        plot([0 0], -[0.05 0.25], 'k', 'LineWidth', S.lwMain)
        plot([alpha alpha], -[0.05 0.25], 'k', 'LineWidth', S.lwMain)
        x_ruler = linspace(0,alpha,100);
        y_ruler = -0.15*ones(size(x_ruler));
        plot(x_ruler, y_ruler, 'k', 'LineWidth', S.lwMain)
        text(alpha/2+0.15, -1.2, ...
            '$\alpha$', ...
            'Interpreter','latex', 'FontSize',S.fsAnnot, ...
            'HorizontalAlignment','center', 'VerticalAlignment','top');
        plot([alpha/2 alpha/2+0.15],[-0.15 -1.3],'k','LineWidth',S.lwLead)

        % dimension ruler for U1
        plot([0 0], -[0.05 0.25], 'k', 'LineWidth', S.lwMain)
        plot([-rho -rho], -[0.05 0.25], 'k', 'LineWidth', S.lwMain)
        x_ruler = linspace(-rho,0,50);
        y_ruler = -0.15*ones(size(x_ruler));
        plot(x_ruler, y_ruler, 'k', 'LineWidth', S.lwMain)
        text(-rho/2-0.15, -1.2, ...
            '$\rho$', ...
            'Interpreter','latex', 'FontSize',S.fsAnnot, ...
            'HorizontalAlignment','center', 'VerticalAlignment','top');
        plot([-rho/2 -rho/2-0.15],[-0.15 -1.25],'k','LineWidth',S.lwLead)
    end

    %% --- a*(t) label ---
    label_positions = [
                        [1.5 2]
                        [alpha/2+0.2 -1.3]
                        [alpha/2+0.2 -1.3]
                        [alpha/2-0.2  1.3]
                      ];
    vrt_align = [
        "bottom"
        "top"
        "top"
        "bottom"
    ];
    label_position = label_positions(k,:);
    text(label_position(1), label_position(2), ...
        '${\rm a}^{\star}(t)/r$', ...
        'Interpreter','latex', 'FontSize',S.fsAnnot, ...
        'HorizontalAlignment','center', 'VerticalAlignment',vrt_align(k));
    plot([alpha/2 label_position(1)],[0 label_position(2)],'k','LineWidth',S.lwLead)

    %% --- Label region D(alpha), outside the disk with a leader to its edge ---
    pD0 = [alpha + 0.85*cos(D_phi(k)), 0.85*sin(D_phi(k))];     % point just inside the disk boundary
    pD1 = D_pos(k,:);                                           % top-centre of the label
    plot([pD0(1) pD1(1)], [pD0(2) pD1(2)+0.1], 'k', 'LineWidth', S.lwLead);  % stop short of the text
    text(pD1(1), pD1(2), '$\mathcal{D}(\alpha)$', ...
        'Interpreter','latex', ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'FontSize',S.fsAnnot);

    hold off

    %% --- Limits, r-scaled ticks, shared styling, size ---
    ax = gca;
    axis(ax,'equal');

    L = 4;
    xlim(ax,[-L L]);
    ylim(ax,[-L L]);

    ticks = (-ceil(L):ceil(L));
    xticks(ax, ticks);
    yticks(ax, ticks);
    xticklabels(ax, arrayfun(@(x) sprintf('%g', x), ticks, 'UniformOutput', false));
    yticklabels(ax, arrayfun(@(x) sprintf('%g', x), ticks, 'UniformOutput', false));

    tacAxes(ax, S);
    xtickangle(0)

    fig.Units = 'inches';
    fig.Position = [1 1 S.panelW S.panelW];   % one cell of the 2x2 grid (square)

    %% --- Export ---
    if save_figure
        figure_name = ['figure_2_' letter(k) '.pdf'];
        exportgraphics(fig, fullfile('figures',figure_name), 'ContentType','vector');
    end
end

end
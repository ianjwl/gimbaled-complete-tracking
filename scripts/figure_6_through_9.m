function figure_6_through_9(sim, save_figures)
%% ======================================================================
%  figure_6_through_9.m  --  the four time-history figures of Section VI
%
%  Produces four separate figure windows, exported as four PDFs:
%
%      figure_6   actual and reference position, one curve per East,
%                 North, Up component, in a 2:3 axes box
%      figure_7   position-error norm, ||x(t) - x*(t)||
%      figure_8   normalized thrust Tbar(t) = lambda(t) u1(t)/m, one curve
%                 per component, normalized by vehicle weight
%      figure_9   angular offset theta(t) against the bound
%                 theta*(||a*(t)||/r)
%
%  Every panel shades the commanded free-fall window of the ballistic arc
%  in gray and labels it. Figures 6, 7, 8, and 9 share the same time axis
%  limits and the same axis styling so they stack cleanly in the paper.
%
%  Encoding convention: actual signals are colored solids, with the color
%  coding the component; reference signals are black dashes carrying a
%  single legend entry. This keeps the target visible where it coincides
%  with the actual signal, which is nearly everywhere except the agile
%  maneuvers.
%
%  Inputs:
%      sim            struct of simulation output from RUN_ALL; uses
%                         sim.t                         time          [Nx1]
%                         sim.ENU_positions             actual pos.   [Nx3]
%                         sim.ENU_reference_positions   reference     [Nx3]
%                         sim.ENU_Thrust                Tbar(t)       [Nx3]
%                         sim.ENU_Thrust_star           Tbar*(t)      [Nx3]
%                         sim.theta                     theta(t)      [Nx1]
%                         sim.theta_star                theta*(t)     [Nx1]
%      save_figures   (optional, default false) logical; when true, writes
%                     figure_6.pdf through figure_9.pdf
%
%  Notes:
%      The thrust is normalized by W = m g with m = 2.09 kg hard-coded
%      here; if the mass changes in the parameters block, this constant
%      must change with it. Reading it from the parameters struct would
%      remove the coupling.
%
%      The free-fall band is the interval [14.5, 17] s, hard-coded to the
%      reference history of the paper.
%
%  Called by RUN_ALL.
%% ======================================================================

if nargin < 2
    save_figures = false;
end

%% Load workspace...
t = sim.t;
ENU_positions = sim.ENU_positions;
ENU_reference_positions = sim.ENU_reference_positions;
ENU_Thrust = sim.ENU_Thrust;
ENU_Thrust_star = sim.ENU_Thrust_star;
theta = sim.theta;
theta_star = sim.theta_star;

%% ---- Global style: LaTeX everywhere -> Computer Modern -----------------
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% ---- Palette (Okabe–Ito): one color per ENU component -----------------
cE = [0   114 178]/255;   % blue         -> East
cN = [213  94   0]/255;   % vermillion   -> North
cU = [0   158 115]/255;   % bluish green -> Up
compColor = [cE; cN; cU];

colW  = 3.5;    % IEEE single-column width [in]
lwa   = 1.5;    % actual line width
lw    = 1;      % target (black dashed) line width
axLW  = 0.75;   % axis line width
fsAx  = 8;      % tick labels
fsLab = 9;      % axis (x) label
fsLeg = 8;      % legend
fsBand= 7;      % "commanded free-fall" text
fsUnit= 7;      % units text

% Commanded free-fall highlight window
ff      = [14.5 17];
ffColor = [0.70 0.74 0.80];
ffAlpha = 0.35;
tL      = [t(1) t(end)];

% Parameters for thrust normilization by weight
m = sim.m; % kg
g = sim.g; % m/s^2
W = m*g;

%% ---- Legend labels ----------------------------------------------------
lblX = {'${\rm x}_{\rm E}$','${\rm x}_{\rm N}$','${\rm x}_{\rm U}$'};
lblT = {'$\bar{T}_{\rm E}$','$\bar{T}_{\rm N}$','$\bar{T}_{\rm U}$'};

%% ======================================================================
%  Figure 6 : Position history (actual colored solid, target black dashed)
%% ======================================================================
f6 = figure('Theme','light','Units','inches','Position',[1 4 colW 2.5]);
ax6 = axes(f6); hold(ax6,'on');
for k = 1:3
    plot(ax6,t,ENU_positions(:,k),'-','Color',compColor(k,:),'LineWidth',lwa,'DisplayName',lblX{k});
end
plotTargets(ax6,t,ENU_reference_positions,lw,'${\rm x}^{\star}$');
hold(ax6,'off');
xlabel(ax6,'Time [sec]','FontSize',fsLab);
xlim(ax6,tL);
applyTACStyle(ax6,fsAx,axLW);
legend(ax6,'Position',[0.6712 0.6725 0.2199 0.1900],'NumColumns',1,'FontSize',fsLeg);
addFreefall(ax6,ff,ffColor,ffAlpha,'commanded free-fall','above',fsBand);
addUnits(ax6,'m',fsUnit)

%% ======================================================================
%  Figure 7 : Position-error norm  ||x - x*||  (label as NW legend entry)
%% ======================================================================
f7 = figure('Theme','light','Units','inches','Position',[1 1 colW 1.1]);
ax7 = axes(f7); ax7.Units = 'normalized'; ax7.Position = [0.135 0.145 0.835 0.72];
hold(ax7,'on');
enorm = vecnorm(ENU_positions - ENU_reference_positions, 2, 2);
plot(ax7,t,enorm,'-','Color',[0 0 0],'LineWidth',lw,'DisplayName','$\|{\rm x}-{\rm x}^{\star}\|$');
hold(ax7,'off');
xlabel(ax7,'Time [sec]','FontSize',fsLab);
xlim(ax7,tL);
applyTACStyle(ax7,fsAx,axLW);
legend(ax7,'Position',[0.2267 0.8003 0.2329 0.16],'FontSize',fsLeg,'IconColumnWidth',20);
addFreefall(ax7,ff,ffColor,ffAlpha,'commanded free-fall','above',fsBand);
addUnits(ax7,'m',fsUnit)

%% ======================================================================
%  Figure 8 : Target vs. true thrust  (T = lambda u1)
%% ======================================================================
f8 = figure('Theme','light','Units','inches','Position',[5 4 colW 2.5]);
ax8 = axes(f8); hold(ax8,'on');
for k = 1:3
    plot(ax8,t,ENU_Thrust(:,k)/W,'-','Color',compColor(k,:),'LineWidth',lwa,'DisplayName',lblT{k});
end
plotTargets(ax8,t,ENU_Thrust_star/W,lw,'$\bar{T}^{\star}$');
hold(ax8,'off');
xlabel(ax8,'Time [sec]','FontSize',fsLab);
xlim(ax8,tL);
applyTACStyle(ax8,fsAx,axLW);
legend(ax8,'Position',[0.1426 0.645 0.1948 0.2000],'NumColumns',1,'FontSize',fsLeg);
addFreefall(ax8,ff,ffColor,ffAlpha,'commanded free-fall','above',fsBand);
addUnits(ax8,'g',fsUnit)

%% ======================================================================
%  Figure 9 : Angular offset (pi/4 ticks, ylim [0 2*pi], both solid)
%% ======================================================================
f9 = figure('Theme','light','Units','inches','Position',[1 1 colW 1]);
ax9 = axes(f9); ax9.Units = 'normalized'; ax9.Position = [0.135 0.145 0.835 0.72];
hold(ax9,'on');
plot(ax9,t,theta_star,'-','Color',[0 0 0],'LineWidth',lw,'DisplayName','$\theta^{\star}$');
plot(ax9,t,theta,     '-','Color',cE,'LineWidth',lw,'DisplayName','$\theta$');
hold(ax9,'off');
xlabel(ax9,'Time [sec]','FontSize',fsLab);
xlim(ax9,tL);
ylim(ax9,[0 2*pi]);
yticks(ax9,0:pi/2:2*pi);
yticklabels(ax9,{'$0$','$\pi/2$','$\pi$', ...
                 '$3\pi/2$','$2\pi$'});
applyTACStyle(ax9,fsAx,axLW);
legend(ax9,'Position',[0.1447 0.5538 0.1954 0.2049],'FontSize',fsLeg);
addFreefall(ax9,ff,ffColor,ffAlpha,'commanded free-fall','above',fsBand);
addUnits(ax9,'rad',fsUnit)

%% ---- Vector export -----------------------------------------------------
if save_figures
    exportgraphics(f6,'figures\figure_6.pdf','ContentType','vector');
    exportgraphics(f7,'figures\figure_7.pdf','ContentType','vector');
    exportgraphics(f8,'figures\figure_8.pdf','ContentType','vector');
    exportgraphics(f9,'figures\figure_9.pdf','ContentType','vector');
end

end

%% ======================================================================
%  Local functions
%% ======================================================================
function plotTargets(ax,t,M,lwT,label)
% Black dashed target for every ENU component; only the first carries a
% legend entry (they all share one style = "target").
    for k = 1:size(M,2)
        if k == 1
            plot(ax,t,M(:,k),'--','Color','k','LineWidth',lwT,'DisplayName',label);
        else
            plot(ax,t,M(:,k),'--','Color','k','LineWidth',lwT,'HandleVisibility','off');
        end
    end
end

function applyTACStyle(ax,fsAx,axLW)
    grid(ax,'on');
    ax.GridAlpha  = 0.15;
    ax.LineWidth  = axLW;
    ax.FontSize   = fsAx;
    ax.Box        = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.XMinorGrid = 'off';
end

function h = addFreefall(ax,xr,faceColor,faceAlpha,label,mode,fs)
    if nargin < 6, mode = 'inside'; end
    if nargin < 7, fs = 7; end
    yl = ylim(ax); ylim(ax,yl);
    h = patch(ax,[xr(1) xr(2) xr(2) xr(1)],[yl(1) yl(1) yl(2) yl(2)], ...
        faceColor,'FaceAlpha',faceAlpha,'EdgeColor','none','HandleVisibility','off');
    uistack(h,'bottom');
    if ~isempty(label)
        switch mode
            case 'above'
                text(ax,mean(xr),yl(2),label,'Interpreter','latex','FontSize',fs, ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom','Clipping','off');
            otherwise
                text(ax,mean(xr),yl(2),label,'Interpreter','latex','FontSize',fs, ...
                    'HorizontalAlignment','center','VerticalAlignment','top', ...
                    'BackgroundColor','w','Margin',0.5,'EdgeColor',[0.75 0.75 0.75]);
        end
    end
end

function addUnits(ax,label,fs)
    if nargin < 3, fs = 7; end
    xl = xlim(ax); xlim(ax,xl);
    yl = ylim(ax); ylim(ax,yl);
    text(ax,xl(1),yl(2),label,'Interpreter','latex','FontSize',fs, ...
                    'HorizontalAlignment','left','VerticalAlignment','bottom','Clipping','off');
end
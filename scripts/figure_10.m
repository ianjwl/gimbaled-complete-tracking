function figure_10(sim, save_figures)
%% ======================================================================
%  figure_10.m  --  first flip rendered against the uprightness history
%
%  Produces one figure window with two stacked tiles:
%
%      top      the first flip, rendered as a sequence of multirotor poses
%               over the flip window with a thin time-connector line from
%               each pose down to a time axis
%      bottom   the uprightness (-e3)'*lambda(t), with the inverted and
%               upright excursions labeled and the free-fall band shaded
%
%  The two tiles share a time axis so the reader can match each rendered
%  pose to the moment in the uprightness trace. Dashed connector lines are
%  drawn in an overlay axes behind both tiles.
%
%  The poses are drawn by plotObjectTrajectory with the mesh quad.stl.
%  Both belong to third_party/FlightTest-Tools-main and are not part of
%  this repository's own code; RUN_ALL puts that folder on the path, so
%  this function should NOT call addpath itself.
%
%  Inputs:
%      sim            struct of simulation output from RUN_ALL; uses
%                         sim.t                    time history          [Nx1]
%                         sim.attitudes            attitude quaternions   [Nx4]
%                         sim.uprightness          (-e3)'*lambda(t)       [Nx1]
%      save_figures   (optional, default false) logical; when true, writes
%                     figure_10.pdf
%
%  Notes:
%      The flip window [7.85, 8.85] s and the free-fall band [14.5, 17] s
%      are hard-coded to the reference history of the paper, as are the
%      pose instants and the pose-spacing offsets that keep the rendered
%      meshes from overlapping.
%
%      The NED position is displayed in an East-North-Up frame, as in
%      figure_5.
%
%      Figure styling is currently set locally, with local applyTACStyle,
%      addFreefall, and addFlip helpers duplicated from figure_6_through_9.
%
%  Requires on the MATLAB path:
%      plotObjectTrajectory and quad.stl (third_party).
%
%  Called by RUN_ALL.
%% ======================================================================

if nargin < 2
    save_figures = false;
end

%% Initialize plotting variables...
t = sim.t;
attitude_quaternion = sim.attitudes;
uprightness = sim.uprightness;

%% Initialize drawing parameters...
stl_times = [8 linspace(8.06,8.10,4) 8.375 linspace(8.565,8.60,4) 8.75];
stl_positions = linspace(min(stl_times),max(stl_times),numel(stl_times));
stl_positions(2) = stl_positions(2) - 0.01;
stl_positions(7) = stl_positions(7) - 0.01;
stl_inds = nan(size(stl_times));
for ind = 1:numel(stl_times)
    stl_inds(ind) = find(t>=stl_times(ind),1,'first');
end

%% ---- Global style: LaTeX everywhere -> Computer Modern -----------------
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% ---- Palette (Okabe–Ito): one color per ENU component -----------------
c1 = [0   114 178]/255;   % blue         -> East
c2 = [213  94   0]/255;   % vermillion   -> North
c3 = [0   158 115]/255;   % bluish green -> Up
compColor = [c1; c2; c3];

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
ffColor = [0.80 0.84 0.90];
ffAlpha = 0.35;
tL      = [t(1) t(end)];

% Flip highlight window
fl      = [7.85 8.85];
flColor = [0.70 0.74 0.80];
flAlpha = 0.35;

%% ======================================================================
%  Figure Setup
%% ======================================================================
fig = figure('Theme','light','Units','inches','Position',[1 2 colW 2.5]);
tiledlayout(2,1,"TileSpacing","compact")

%% Plot flip
nexttile(1)
scaleFactor = 0.1;
numObjects = numel(stl_inds);
stl = 'quad.stl';
translations = [zeros(size(stl_positions(:))) stl_positions(:) zeros(size(stl_positions(:)))];
rotations = attitude_quaternion(stl_inds,:);
plotObjectTrajectory(translations,rotations,scaleFactor,numObjects,stl,1:numel(stl_inds))
delete(fig.Children.Children.Children(end));
xlabel('Time [sec]')
ylabel('')
zlabel('')
view([0 0])
xlim([7.85 8.85])
lim = 0.065;
ylim([-lim lim])
zlim([-lim lim])
xticks(7.9:0.1:8.8)
yticks([])
zticks([])
yticklabels([])
zticklabels([])
grid on
axis equal
box on
for i = numel(fig.Children.Children.Children):-1:1
    if strcmp(fig.Children.Children.Children(i).Type,'light')
        delete(fig.Children.Children.Children(i))
    end
end
for i = 1:numel(fig.Children.Children.Children)
    delete(fig.Children.Children.Children(i).Children.Children.Children(1))
    delete(fig.Children.Children.Children(i).Children.Children.Children(1))
    delete(fig.Children.Children.Children(i).Children.Children.Children(1))
end
lightangle(35,20)
hold on
for i = 1:numel(stl_inds)
    plot3([stl_positions(i) stl_times(i)],[0 0],[0 -lim],'Color',c1,'LineWidth',1)
end
hold off
applyTACStyle(gca,fsAx,axLW);
% set(gca,'Visible','off')
% set(gca,'XColor', 'none','YColor','none')

%% Plot uprightness
nexttile(2)
plot(t,uprightness,'-','Color',[0 0 0],'LineWidth',lw,'DisplayName','$(-{\rm e}_3)^\top\lambda(t)$');
xlabel('Time [sec]','FontSize',fsLab);
xlim(tL);
ylim([-1.1 1.1])
yticks(-1:1:1)
applyTACStyle(gca,fsAx,axLW);
legend('Position',[0.6712 0.1952 0.26 0.0810],'FontSize',fsLeg,'IconColumnWidth',20);
addFreefall(gca,ff,ffColor,ffAlpha,'commanded free-fall','above',fsBand);
text(tL(1)+0.2,-1.0,'Inverted','Interpreter','latex','FontSize',fsBand, ...
                    'HorizontalAlignment','left','VerticalAlignment','bottom','Clipping','off');
text(tL(1)+0.2, 1.0,'Upright', 'Interpreter','latex','FontSize',fsBand, ...
                    'HorizontalAlignment','left','VerticalAlignment','top','Clipping','off');
addFlip(gca,fl,flColor,flAlpha,'first flip','above',fsBand);

%% Add annotation lines
% a1 = annotation(fig,"line",[0.13 0.375],[0.675 0.51],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',1);
% uistack(a1, 'bottom');
% a2 = annotation(fig,"line",[0.9 0.4],[0.675 0.51],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',1);
% uistack(a2, 'bottom');

%% Add connector lines in a background axes (behind the tiles, incl. tick labels)
bg = axes(fig,'Units','normalized','Position',[0 0 1 1], ...
          'XLim',[0 1],'YLim',[0 1],'Color','none', ...
          'HitTest','off','HandleVisibility','off');
hold(bg,'on');
plot(bg,[0.13 0.375],[0.675 0.51],'--','Color',[0.85 0.85 0.85],'LineWidth',1);
plot(bg,[0.90 0.400],[0.675 0.51],'--','Color',[0.85 0.85 0.85],'LineWidth',1);
axis(bg,'off');
uistack(bg,'bottom');   % send the whole overlay behind the tiled axes

%% Export
if save_figures
    exportgraphics(fig,'figures\figure_10.pdf', 'ContentType','vector');
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

function h = addFlip(ax,xr,faceColor,faceAlpha,label,mode,fs)
    if nargin < 6, mode = 'inside'; end
    if nargin < 7, fs = 7; end
    yl = ylim(ax); ylim(ax,yl);
    h = patch(ax,[xr(1) xr(2) xr(2) xr(1)],[yl(1) yl(1) yl(2) yl(2)], ...
        faceColor,'FaceAlpha',faceAlpha,'EdgeColor','none','HandleVisibility','off');
    uistack(h,'bottom');
    if ~isempty(label)
        switch mode
            case 'above'
                text(ax,min(xr),mean(yl),label,'Interpreter','latex','FontSize',fs, ...
                    'HorizontalAlignment','right','VerticalAlignment','bottom','Clipping','off');
            otherwise
                text(ax,min(xr),mean(yl),label,'Interpreter','latex','FontSize',fs, ...
                    'HorizontalAlignment','right','VerticalAlignment','top', ...
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
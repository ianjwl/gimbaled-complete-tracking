function figure_5(sim, save_figures)
%% ======================================================================
%  figure_5.m  --  three-dimensional trajectory with rendered vehicle poses
%
%  Draws the flown trajectory through the obstacle course with the engine
%  reference overlaid, and renders the multirotor mesh at a set of instants
%  so the two flips and the ballistic arc are visible in the trajectory
%  itself. Also draws the flight environment: the takeoff polygon, the
%  floor, and the two obstacle hoops.
%
%  Signals: actual position x(t) as a colored solid curve, reference
%  position x*(t) as a black dashed curve. The dashes are black rather than
%  color-matched so the reference stays legible where it lies on top of the
%  actual, and so both survive a grayscale print.
%
%  The poses are drawn by plotObjectTrajectory with the mesh quad.stl.
%  Both belong to third_party/FlightTest-Tools-main and are not part of
%  this repository's own code; RUN_ALL puts that folder on the path, so
%  this function should NOT call addpath itself.
%
%  The pose instants and the dimensions of the hoops, floor, and takeoff
%  polygon are hard-coded to the reference history of the paper. If the
%  reference history is changed, the environment drawn here no longer
%  corresponds to the trajectory and this function must be edited too.
%
%  Inputs:
%      sim            struct of simulation output from RUN_ALL; uses
%                         sim.t                    time history          [Nx1]
%                         sim.positions            inertial pos., NED     [Nx3]
%                         sim.reference_positions  reference pos., NED    [Nx3]
%                         sim.attitudes            attitude quaternions   [Nx4]
%      save_figures   (optional, default false) logical; when true, writes
%                     figure_5.pdf
%
%  Notes:
%      The NED position is displayed in an East-North-Up frame by plotting
%      [North, East, -Down] against the axes labeled Easting, Northing,
%      Altitude. The same component swap appears in ROOT's postprocessing.
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
position = sim.positions;
desired_position = sim.reference_positions;
attitude_quaternion = sim.attitudes;
t = sim.t;

%% Initialize drawing parameters...
stl_times = [0.5 1 2 4 6 7.2 9.35 10.7 14 14.25 14.5 14.75 15 15.25 15.5 15.75 16 16.25 16.5 16.75 17 17.25 17.5 18 20 21.5 22.5 23.5];
stl_inds = nan(size(stl_times));
for ind = 1:numel(stl_times)
    stl_inds(ind) = find(t>=stl_times(ind),1,'first');
end
text_h_alignments = [repmat("left",1,6) "center" "left" repmat("r",1,5) "left" "left" "center" repmat("r",1,12) "right"];
text_v_alignments = ["middle" "middle" "top" "top" "middle" "middle" "bottom" repmat("middle",1,8) "bottom" repmat("middle",1,8) "top" "top" "middle" "middle" "middle"];
text_h_offsets = [repmat(0.75,1,4) 0.25 0.75 0 0.75 -0.25 -0.75 -0.5 -0.75 -0.75 0.5 0.5 0 -0.625 -0.625 -0.5 -0.875 -0.45 -0.75 -1 -1 -0.625 -0.625 -0.75 -0.75 -0.25];
text_v_offsets = [0 0 -0.25 -0.25 -0.25 0 1 0 -0.25 repmat(0,1,6) 1 repmat(0,1,8) -0.25 -0.25 -0.125 0 0.25];

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

%% ======================================================================
%  Figure Setup
%% ======================================================================
fig = figure('Theme','light','Units','inches','Position',[1 2 colW 7]);
% set(gcf,'Color','white')
scaleFactor = 1;
numObjects = 20;
stl = 'quad.stl';
translations = position;
rotations = attitude_quaternion;
plotObjectTrajectory(translations,rotations,scaleFactor,numObjects,stl,stl_inds(1:end-3))
delete(fig.Children.Children(end));
grid on
axis equal
box on
xlabel('\textbf{Easting [m]}','Interpreter','latex','FontSize',fsLab)%,'Position',[120 -120]
ylabel('\textbf{Northing [m]}','Interpreter','latex','FontSize',fsLab)%,'Position',[310 5]
text(-3.8,-2,19,'\textbf{Altitude [m]}','Interpreter','latex','FontSize',fsLab)
% zlabel('\textbf{Altitude [m]}','Interpreter','latex')
% xlabel('')
% ylabel('')
zlabel('')
xlim([-4 4])
ylim([-2 6])
zlim([0 20])
xticks(-4:2:4)
yticks(-2:2:6)
zticks(0:2:20)
view([21.5 10]) % view([35 9])
fig.Children.XLabel.Position = [0.3 -2 -0.875];
fig.Children.YLabel.Position = [5 1.5 -1.2];
hold on
plot3(position(:,2),position(:,1),-position(:,3),'LineWidth',lwa,'Color',c3,'DisplayName','${\rm x}(t)$')
plot3(desired_position(:,2),desired_position(:,1),-desired_position(:,3),'--k','LineWidth',lw,'DisplayName','${\rm x}^{\star}(t)$')
hold off
for i = numel(fig.Children.Children):-1:1
    if strcmp(fig.Children.Children(i).Type,'light')
        delete(fig.Children.Children(i))
    end
end
lightangle(35,60)
stl_times(1) = 0; % I'm plotting 0.5s since 0s looks like it's in the floor, but the label should reflect 0s
stl_inds = [stl_inds stl_inds(5)]; % Adding in a second time callout at t = 13
stl_times = [stl_times 13]; 
for i = 1:numel(stl_inds)
    ti = stl_times(i);
    ii = stl_inds(i);
    Ei = position(ii,2) + text_h_offsets(i);
    Ni = position(ii,1);
    Ui =-position(ii,3) + text_v_offsets(i);
    text(Ei,Ni,Ui,ti+"s",'FontSize',fsBand,...
        'HorizontalAlignment',text_h_alignments(i),...
        'VerticalAlignment',"middle")%text_v_alignments(i)
end
legend('Location','northwest','FontSize',fsLeg)

hold on
quiver3(2,4,19.25,0.8,0,-1,'off','filled','Color','k','HandleVisibility','off','MaxHeadSize',10)
hold off
text(2.125,4,19.5,'Flips not visible','HorizontalAlignment','right')

%% ======================================================================
%  Environment
%% ======================================================================
hold on
starting_polygon = fill3(0.5*cos(linspace(0,2*pi,100)), 0.5*sin(linspace(0,2*pi,100)), zeros(1,100),c1,'FaceAlpha', 0.5,'HandleVisibility','off');
floor_triangle = plot3([0 3 -3 0],[0 4 4 0],[0 0 0 0],'Color',c1,'HandleVisibility','off');
right_loop = plot3(3+cos(linspace(0,2*pi,100)), 4+sin(linspace(0,2*pi,100)), 3.5*ones(1,100),'Color',c1,'LineWidth',2,'HandleVisibility','off');
right_loop_1st_leg = plot3([4 4],[4 4],[0 3.5],'Color',c1,'LineWidth',2,'HandleVisibility','off');
right_loop_2nd_leg = plot3([2 2],[4 4],[0 3.5],'Color',c1,'LineWidth',2,'HandleVisibility','off');
left_loop = plot3(-3+cos(linspace(0,2*pi,100)), 4+sin(linspace(0,2*pi,100)), 3.5*ones(1,100),'Color',c1,'LineWidth',2,'HandleVisibility','off');
left_loop_1st_leg = plot3([-2 -2],[4 4],[0 3.5],'Color',c1,'LineWidth',2,'HandleVisibility','off');
left_loop_2nd_leg = plot3([-4 -4],[4 4],[0 3.5],'Color',c1,'LineWidth',2,'HandleVisibility','off');
hold off

%% Export
if save_figures
    exportgraphics(fig,'figures\figure_5.pdf', 'ContentType','vector');
end

end
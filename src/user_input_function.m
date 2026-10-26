function user_input = user_input_function(t,vehicle_state,parameters)
%USER_INPUT_FUNCTION  Builds the user command for the control law.
%   user_input = USER_INPUT_FUNCTION(t,vehicle_state,parameters)
%   returns the 6x1 user command consumed by CONTROL_LAW:
%
%       user_input(1:3)   reference position history x*(t)
%       user_input(4:6)   desired pointing vector nu(t)
%
%   The active reference is the obstacle-course history of
%   ALTERNATE_REFERENCE_FUNCTION, which is the history used in the paper's
%   simulation. It already returns x_user in ENU order as [Y; X; -Z], so
%   user_input(1:3) is ENU while the vehicle state is NED; the plotting
%   code in RUN_ALL performs the matching component swap for display.
%
%   The pointing vector is prescribed as
%       nu(t) = [v_N(t); v_E(t); 0]
%   i.e. the ground-track velocity with the vertical component removed.
%   This is the signal that the second-order low-pass filter of
%   eq. (outer_claw_2) acts on before it reaches the inner loop, whose
%   auxiliary-direction objective drives mu(t) to span{a*(t),nu(t)}.
%
%   When parameters.user_evaluates_error is true, the command is overridden
%   if the vehicle has drifted more than
%   parameters.user_position_error_threshold from the reference, in which
%   case the reference is replaced by a point 0.8 times the current
%   position and nu is set to e_1. This is a debugging aid; RUN_ALL sets
%   user_evaluates_error to false.
%
%   The eight commented-out blocks above the active branch are alternative
%   reference histories (L path, wiggle, vertical sinusoid, planar circle,
%   vertical helix, ladel, and two obstacle-course variants) retained for
%   reference. They are not executed. Several call SSAT and the last of
%   them uses the not_free_falling output of ALTERNATE_REFERENCE_FUNCTION.
%
%   Inputs:
%       t              current time, scalar
%       vehicle_state  vehicle state, layout as in MULTIROTOR_FLOW_MAP;
%                      uses vehicle_state(1:3) only under the error override
%                      and vehicle_state(4:6) for nu
%       parameters     struct of system and controller parameters; uses
%                      parameters.user_evaluates_error and
%                      parameters.user_position_error_threshold
%
%   Outputs:
%       user_input     6x1 command, layout above
%
%   Requires on the MATLAB path:
%       alternate_reference_function (active branch).
%
%   Called by CONTROL_LAW_FOR_SIMULATION, MULTIROTOR_JUMP_GUARD, and
%   MULTIROTOR_LYAPUNOV_FUNCTION.

%% Obstacle course, ground-track pointing (used in the paper)
user_position = alternate_reference_function(t);
user_direction = [vehicle_state(4:5); 0]; 
user_input = [user_position; user_direction];

%% L path
% user_position = [0;0;-10]*max([0 min([(t-0)/1.6 1])]) + [0;0;10]*max([0 min([(t-10)/1.6 1])]) + [10;10;0]*max([0 min([(t-20)/5 1])]);
% user_direction = [vehicle_state(4:5); 0]; % Vehicle commanded to point along NE part of inertial NED velocity
% user_input = [user_position; user_direction];

%% Wiggle path
% user_position = ssat(t,0,5)*([5;0;0]*cos(2*pi*t/10) + [0;5;0]*sin(2*pi*t/10) + [0;0;-1]*sin(2*pi*t/2.5));
% user_direction = [vehicle_state(4:5); 0]; % Vehicle commanded to point along NE part of inertial NED velocity
% user_input = [user_position; user_direction];

%% Vertical Sinusoid
% A = 1.2*9.81; % peak acceleration
% P = 4*5; % period
% H = A/((2*pi/P)^2/2); % height
% user_position = [0;0;-1]*( H*(1-cos(2*pi*t/P))/2 );
% user_direction = [1;0;0];
% user_input = [user_position; user_direction];

%% Planar circle
% P = 10; % Period
% R = 10; % Radius
% user_position = ssat(t,0,1.5*P)*R*( cos(2*pi*t/(P))*[1;0;0] + sin(2*pi*t/(P))*[0;1;0]);
% user_direction = [1;0;0];
% user_input = [user_position; user_direction];

%% Vertical Helix
% A = 1.2*9.81; % Peak Acceleration
% P = 4*5; % Period
% H = A/((2*pi/P)^2/2); % Height
% R = 10; % Radius
% user_position = [0;0;-1]*( H*(1-cos(2*pi*t/P))/2 ) ...
%      + ssat(t,0,1.5*P/4)*R*( cos(2*pi*t/(P/4))*[1;0;0] + sin(2*pi*t/(P/4))*[0;1;0]);
% user_direction = [1;0;0];
% user_input = [user_position; user_direction];

%% Ladel path
% user_position = [0;0;-10]*max([0 min([(t-0)/1.6 1])]) + [0;0;10]*max([0 min([(t-6)/1.6 1])]) ...
%     + 2.5/sqrt(2)*([1;1;0]*(1-cos(2*pi*max([0 min([(t-11)/3 1])]))) + [1;-1;0]*sin(2*pi*max([0 min([(t-11)/3 1])])));
% user_direction = (1-ssat(t,14,15))*[vehicle_state(4:5); 0]/(norm([vehicle_state(4:5); 0])+0.01)*(norm([vehicle_state(4:5); 0])>1) + ssat(t,14,15)*[1;0;0]; 
% user_input = [user_position; user_direction];

%% Obstacle Course
% [x_user,~,~,~,~,not_free_falling] = alternate_reference_function(t);user_position = x_user;
% raw_user_direction = (1-ssat(t,21.5,22))*[vehicle_state(4:5); 0]/(norm([vehicle_state(4:5); 0])+0.01)*(norm([vehicle_state(4:5); 0])>1) + ssat(t,21.5,22)*[1;0;0]; 
% user_direction = not_free_falling*raw_user_direction;
% user_input = [user_position; user_direction];

%% User error evaluation
% if parameters.user_evaluates_error
%     vehicle_position = vehicle_state(1:3);
%     if norm(vehicle_position - user_position) > parameters.user_position_error_threshold
%         user_input = [ 0.8*vehicle_position ; [1;0;0] ];
%     end
% end


























%% Old stuff
% % Assign user input...
% % user_position = [0;0;1.3*9.81/2*(t-5)^2]*(t>5)*(t<10) + [0;0;1.3*9.81/2*(5)^2 - 1.3*9.81/2*(t-10)^2]*(t>=10)*(t<15) + [5;5;0]*(t>20);
% % user_position = [0;0;-10]*(t>=5)*(t<15) + [5;5;0]*(t>=20);
% % user_position = [0;0;-10]*(t>=5)*(t<15) + [5;5;0]*max([0 min([(t-20)/2.5 1])]);
% % user_position = [0;0;-max([0 min([2*9.81/2*(t-0)*abs(t-0) 10])])] + [0;0;max([0 min([2*9.81/2*(t-10)*abs(t-10) 10])])] + [10;10;0]*max([0 min([(t-20)/5 1])]);
% % user_position = [0;0;-1]*max([0 min([1.2*9.81/2*(t-0)*abs(t-0) 10])]) + [0;0;10]*max([0 min([(t-10)/1.6 1])]) + [10;10;0]*max([0 min([(t-20)/5 1])]);
% user_position = [0;0;-10]*max([0 min([(t-0)/1.6 1])]) + [0;0;10]*max([0 min([(t-10)/1.6 1])]) + [10;10;0]*max([0 min([(t-20)/5 1])]);
% user_direction = [vehicle_state(4:5); 0]; % Vehicle commanded to point along NE part of inertial NED velocity
% user_input = [user_position; user_direction];
% % user_input = 0*user_input;
% 
% % user_position = cos(2*pi*t/10)*[5; 1; 0] + sin(2*pi*t/10)*[-1; 7; 0] + [0; 0; -10];
% % user_direction = [vehicle_state(4:5); 0];
% % user_input = [user_position; user_direction];
% % 
% A = 1.2*9.81; % peak acceleration
% P = 4*5; % period
% H = A/((2*pi/P)^2/2); % height
% 
% % user_position = [0;0;-1]*( H*(1-cos(2*pi*t/P))/2 );
% user_position = ssat(t,0,0.1)*[0;0;-1]*( H*(1-cos(2*pi*t/P))/2 ) ...
%      + ssat(t,40,42)*(H/3)*( cos(2*pi*t/(P/3))*[1;0;0] + sin(2*pi*t/(P/3))*[0;1;0]);
% % user_direction = [vehicle_state(4:5); 0];
% % user_direction = vehicle_state(4:6);
% user_direction = [1;0;0];
% user_input = [user_position; user_direction];
% % 
% % initial_acceleration = 0;
% % maximum_acceleration = 1.2*9.81;
% % inflection_time = 15;
% % 
% % user_position = [0;0;-1]*maximum_acceleration*t.^3*(t - 4*inflection_time)/12/inflection_time^2;
% % user_direction = [1;0;0];
% % user_input = [user_position; user_direction];
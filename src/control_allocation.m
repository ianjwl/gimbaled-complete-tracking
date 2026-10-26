function [f,discrete_CA_controller_state_tpdt,outputs] = control_allocation(u_baseline,discrete_CA_controller_state,parameters,outputs)
%CONTROL_ALLOCATION  Maps the collective thrust and torque command to 
%   individual actuator forces.
%   [f,xd_CA_tpdt,outputs] = CONTROL_ALLOCATION(u_baseline,xd_CA, ...
%       parameters,outputs)
%   distributes the 4x1 net command u_baseline = [thrust; torque] across
%   the four rotors and returns the resulting 4x1 vector of individual
%   motor forces f.
%
%   Mixing matrix. M = parameters.mixing_matrix is built in RUN_ALL and maps
%   motor forces to the net wrench, u = M*f, with rows [thrust; roll;
%   pitch; yaw] and one column per motor. For the arrangement used in RUN_ALL
%   (arm length ell, yaw moment coefficient kt) the rows are
%       [  1     1     1     1  ]
%       [-ell   ell   ell  -ell ]/sqrt(2)
%       [-ell  -ell   ell   ell ]/sqrt(2)
%       [  kt   -kt    kt   -kt ]
%   M is square, so the allocation is the plain inverse f = M\u.
%
%   Active path: RUN_ALL sets apply_QP to false, in which case the allocation
%   is the inverse
%       f = M\u_baseline,      u = u_baseline,
%   and the saturation levels are reset to s = [1;1;1;1]. This is the
%   configuration used for the results in the paper: the vehicle is
%   assumed to have sufficient control authority, so the net thrust and
%   torque requested by the inner loop are reproduced exactly by the four
%   rotors and no actuator limit is active.
%
%   Inactive paths: When apply_QP is true, one of seven quadratic-program
%   allocators is selected by parameters.qp_version. They trade thrust
%   accuracy against torque accuracy subject to the per-motor bounds
%   0 <= f <= parameters.max_individual_motor_thrust, and differ in how
%   they prioritize the components of u. All seven call quadprog and so
%   depend on the Optimization Toolbox; none is used for the results in
%   this paper, and the toolbox is not required when apply_QP is false
%   except that RUN_ALL still builds parameters.qp_options at setup.
%
%   Actuator limits. When parameters.saturate_controls is true, f is
%   clamped to [lb,ub] = [0, max_individual_motor_thrust] and u is
%   recomputed as M*f, so after clamping the net command is no longer
%   equal to u_baseline. RUN_ALL sets saturate_controls to false.
%
%   Naming caution: Inside the inner loop, u1 denotes the collective thrust
%   and u2 the angular acceleration. Here u1 and u2 are parsed from
%   u_baseline, so u1 is the collective thrust and u2 is the torque: the 
%   inner loop has already converted its angular acceleration command to a
%   torque. The quantities are related but not identical, and the local
%   names are reused for convenience.
%
%   Controller state. xd_CA is the 4x1 vector of saturation levels s_lp.
%   Here, "lp" means "last pass". It exists to warm-start the QP allocators; 
%   in the active path it is reset to [1;1;1;1] on every call, so it carries 
%   no information from one step to the next. It is kept in the controller 
%   state layout regardless so that the state indexing in CONTROL_LAW is 
%   unchanged if the QP allocators are ever enabled.
%
%   Inputs:
%       u_baseline   4x1 net command [thrust; torque] from the inner loop,
%                    the signal [u1;u2] of the paper
%       xd_CA        4x1 saturation levels s_lp (i.e., s_last_pass)
%       parameters   struct of system and controller parameters; uses
%                    parameters.mixing_matrix,
%                    parameters.max_individual_motor_thrust, and
%                    parameters.saturate_controls, plus, only when
%                    apply_QP is true, parameters.apply_QP,
%                    parameters.qp_version, parameters.qp_slack,
%                    parameters.lambda_qp, and parameters.qp_options
%       outputs      diagnostics struct, threaded through
%
%   Outputs:
%       f            4x1 individual motor forces
%       xd_CA_tpdt   4x1 saturation levels at t + dt
%       outputs      diagnostics struct, extended with outputs.f_baseline
%                    (the inverse allocation of the input command),
%                    outputs.f, and outputs.u (the net command actually
%                    realized after any clamping)
%
%   Requires on the MATLAB path:
%       unit_vector, and quadprog only when apply_QP is true.
%
%   Called by CONTROL_LAW.

% Load parameter values...
mixing_matrix = parameters.mixing_matrix;
max_individual_motor_thrust = parameters.max_individual_motor_thrust;
apply_QP = parameters.apply_QP;
qp_slack = parameters.qp_slack;
lambda_qp = parameters.lambda_qp;
qp_version = parameters.qp_version;
saturate_controls = parameters.saturate_controls;

% Parse state...
s_lp = discrete_CA_controller_state';

% Rename saturation parameters...
M = mixing_matrix;
fmin = 0;
fmax = max_individual_motor_thrust;
lb = fmin*ones(4,1);
ub = fmax*ones(4,1);

% Assemble control data...
f_baseline = (M\u_baseline);
f = f_baseline;
u = u_baseline;
u1 = u(1,1);
u2 = u(2:4,1);

% Manage QP solution variables...
s1 = s_lp(1);
s2 = s_lp(2);
s3 = s_lp(3);
s4 = s_lp(4);

% Solve control allocation QPs...
if apply_QP
    if qp_version == 1
        
        % First QP...
        Hqp1 = diag([0 1]);
        fqp1 = -[0 1]';
        temp = (M\[u1 0; zeros(3,1) u2]);
        Aqp1 = [temp; -temp];
        bqp1 = [ub; -lb] - qp_slack;
        sqp1 = quadprog(Hqp1,fqp1,Aqp1,bqp1,[],[],[0;0],[],s_lp(1:2),parameters.qp_options);
        s2 = sqp1(2);
    
        % Second QP...
        Hqp2 = 1;
        fqp2 = -1;
        Aqp2 = Aqp1*unit_vector(1,2);
        bqp2 = bqp1 - Aqp1*unit_vector(2,2)*s2 + qp_slack;
        s1 = quadprog(Hqp2,fqp2,Aqp2,bqp2,[],[],0,[],sqp1(1),parameters.qp_options);

        % Assemble control output...
        u = [s1*u1; s2*u2];
        f = (M\u);

        % Set s3 and s4 to avoid errors...
        s3 = 1;
        s4 = 1;
    
    elseif qp_version == 2

        % Solve a single QP...  
        Hqp = diag([lambda_qp(1) 1]);
        fqp = -Hqp*[1;1];
        temp = (M\[u1 0; zeros(3,1) u2]);
        Aqp = [temp; -temp];
        bqp = [ub; -lb];
        sqp = quadprog(Hqp,fqp,Aqp,bqp,[],[],[0;0],[],s_lp(1:2),parameters.qp_options);
        s1 = sqp(1);
        s2 = sqp(2);

        % Assemble control output...
        u = [s1*u1; s2*u2];
        f = (M\u);

        % Set s3 and s4 to avoid errors...
        s3 = 1;
        s4 = 1;

    elseif qp_version == 3

        % Separate control components...
        u_thrust = u1;
        u_roll_and_pitch = u2(1:2);
        u_yaw = u2(3);
        
        % First QP...
        Hqp1 = diag([0 1 0]);
        fqp1 = -[0 1 0]';
        temp = (M\[u_thrust 0 0; zeros(2,1) u_roll_and_pitch zeros(2,1); 0 0 u_yaw]);
        Aqp1 = [temp; -temp];
        bqp1 = [ub; -lb] - 2*qp_slack;
        sqp1 = quadprog(Hqp1,fqp1,Aqp1,bqp1,[],[],[0;0;0],[],s_lp(1:3),parameters.qp_options);
        s2 = sqp1(2);
        
        % Second QP...
        Hqp2 = diag([1 0]);
        fqp2 = -[1 0]';
        Aqp2 = Aqp1(:,[1 3]);
        bqp2 = bqp1 + qp_slack - Aqp1(:,2)*s2;
        sqp2 = quadprog(Hqp2,fqp2,Aqp2,bqp2,[],[],[0;0],[],[sqp1(1); sqp1(3)],parameters.qp_options);
        s1 = sqp2(1);
    
        % Third QP...
        Hqp3 = 1;
        fqp3 = -1;
        Aqp3 = Aqp2(:,2);
        bqp3 = bqp2 + qp_slack - Aqp2(:,1)*s1;
        s3 = quadprog(Hqp3,fqp3,Aqp3,bqp3,[],[],0,[],sqp2(2),parameters.qp_options);

        % Assemble control output...
        u = [s1*u_thrust; s2*u_roll_and_pitch; s3*u_yaw];
        f = (M\u);

        % Set s4 to avoid errors...
        s4 = 1;

    elseif qp_version == 4

        % Separate control components...
        u_thrust = u1;
        u_roll_and_pitch = u2(1:2);
        u_yaw = u2(3);

        % Solve a single QP...  
        Hqp = diag([lambda_qp(1) 1 lambda_qp(2)]);
        fqp = -Hqp*[1;1;1];
        temp = (M\[u_thrust 0 0; zeros(2,1) u_roll_and_pitch zeros(2,1); 0 0 u_yaw]);
        Aqp = [temp; -temp];
        bqp = [ub; -lb];
        sqp = quadprog(Hqp,fqp,Aqp,bqp,[],[],[0;0;0],[],s_lp(1:3),parameters.qp_options);
        s1 = sqp(1);
        s2 = sqp(2);
        s3 = sqp(3);

        % Assemble control output...
        u = [s1*u_thrust; s2*u_roll_and_pitch; s3*u_yaw];
        f = (M\u);

        % Set s4 to avoid errors...
        s4 = 1;

    elseif qp_version == 5
        
        % First QP...
        Hqp1 = M'*diag([0 1 1 1])*M;
        fqp1 = -M'*[zeros(1,3); eye(3)]*u2;
        lbqp1 = lb + qp_slack;
        ubqp1 = ub - qp_slack;
        sqp1 = quadprog(Hqp1,fqp1,[],[],[],[],lbqp1,ubqp1,s_lp,parameters.qp_options);
        u2qp1 = [zeros(3,1) eye(3)]*M*sqp1;
    
        % Second QP...
        Hqp2 = M'*diag([1 0 0 0])*M;
        fqp2 = -M'*[1;0;0;0]*u1;
        Aeqqp2 = [zeros(3,1) eye(3)]*M;
        beqqp2 = u2qp1;
        lbqp2 = lb;
        ubqp2 = ub;
        sqp2 = quadprog(Hqp2,fqp2,[],[],Aeqqp2,beqqp2,lbqp2,ubqp2,sqp1,parameters.qp_options);
        uqp2 = M*sqp2;

        % Assemble control output...
        u = uqp2;
        f = (M\u);

        % Set s1, s2, s3 to avoid errors...
        s1 = sqp2(1);
        s2 = sqp2(2);
        s3 = sqp2(3);
        s4 = sqp2(4);

    elseif qp_version == 6
        
        % First QP...
        Hqp1 = M'*diag([lambda_qp(1) 1 1 1])*M;
        fqp1 = -M'*diag([lambda_qp(1) 1 1 1])*u;
        lbqp1 = lb;
        ubqp1 = ub;
        sqp1 = quadprog(Hqp1,fqp1,[],[],[],[],lbqp1,ubqp1,s_lp,parameters.qp_options);
        uqp1 = M*sqp1;
    
        % Assemble control output...
        u = uqp1;
        f = (M\u);

        % Set s1, s2, s3 to avoid errors...
        s1 = sqp1(1);
        s2 = sqp1(2);
        s3 = sqp1(3);
        s4 = sqp1(4);

    elseif qp_version == 7

        % This QP ONLY enforces lower bounds at 0 on the thrust.  
        lbqp = zeros(4,1);
        ubqp = inf(4,1);
        % This QP enforces the torque command as a constraint and then
        % drives the total thrust to match as closely as possible.  
        
        % First QP...
        Hqp1 = M'*diag([1 0 0 0])*M;
        fqp1 = -M'*[1; zeros(3,1)]*u1;
        Aeqqp1 = [zeros(3,1) eye(3)]*M;
        beqqp1 = u2;
        lbqp1 = lbqp;
        sqp1 = quadprog(Hqp1,fqp1,[],[],Aeqqp1,beqqp1,lbqp1,[],s_lp,parameters.qp_options);
        uqp1 = M*sqp1;

        % Assemble control output...
        u = uqp1;
        f = (M\u);

        % Set s1, s2, s3, s4 to avoid errors...
        s1 = sqp1(1);
        s2 = sqp1(2);
        s3 = sqp1(3);
        s4 = sqp1(4);

    else

        error('control_allocation: QP version not recognized.')

    end
else
    
    % Set s1, s2, s3, and s4...
    s1 = 1;
    s2 = 1;
    s3 = 1;
    s4 = 1;

end

% Apply saturations...
if saturate_controls
    f = (M\u);
    f = max([lb min([f ub],[],2)],[],2);
    u = M*f;
end

% Manage QP solution variables...
s = [s1 s2 s3 s4];

% Calculate s_lp_tpdt...
s_lp_tpdt = s;

% Assemble control-allocation control state increment...
discrete_CA_controller_state_tpdt = s_lp_tpdt(:);

% Manage outputs...
outputs.f_baseline = f_baseline;
outputs.f = f;
outputs.u = u;

end
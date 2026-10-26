function [attitude_matrix,attitude_quaternion_kinematic_matrix] = attitude_quaternion_to_matrix(attitude_quaternion)
%ATTITUDE_QUATERNION_TO_MATRIX  Body-to-inertial rotation matrix from a quaternion.
%   [R,E] = ATTITUDE_QUATERNION_TO_MATRIX(q) returns the rotation matrix R
%   corresponding to the scalar-first unit quaternion q, along with the
%   quaternion kinematic matrix E(q). The two are computed together because
%   the same factorization supplies both.
%
%   Convention. R maps vectors from the body frame to the inertial frame,
%       v_inertial = R*v_body,
%   so the inertial-to-body matrix is R' and is obtained by transposing.
%   The relationship is body-to-inertial, which is the convention used
%   throughout this code; see the remarks in
%   ATTITUDE_QUATERNION_KINEMATIC_MATRIX, which builds E from the same
%   matrices.
%
%   Implementation. Rather than assembling R component by component, the
%   function forms the co-kinematic matrix
%       Eco(q) = (1/2)*[ -qbar'; q0*I - hat(qbar) ]
%   and evaluates
%       R = 4*Eco(q)'*E(q).
%   The two matrices differ only in the sign of the hat(qbar) block, so this
%   is the usual quaternion-to-matrix formula written in factored form, and
%   it reuses the kinematic matrix that the caller often needs anyway.
%   Written out, the equivalent component form is
%       R = (q0^2 - qbar'*qbar)*I + 2*qbar*qbar' + 2*q0*hat(qbar).
%
%   The quaternion is assumed to be a unit quaternion; no normalization is
%   performed here.
%
%   Inputs:
%       attitude_quaternion  4x1 unit quaternion, scalar first
%
%   Outputs:
%       attitude_matrix                    3x3 rotation matrix, body to inertial
%       attitude_quaternion_kinematic_matrix  4x3 matrix E(q), as in
%                                          ATTITUDE_QUATERNION_KINEMATIC_MATRIX
%
%   Requires on the MATLAB path: hat_map.
%
%   Called by MULTIROTOR_FLOW_MAP, MULTIROTOR_JUMP_GUARD,
%   MULTIROTOR_LYAPUNOV_FUNCTION, and INNER_LOOP_CONTROLLER.

%% Quaternion Components
    % This code uses the scalar first convention for qiaternions and
    % abbreviates the scalar part as q0 and the vector part as qbar.  
attitude_quaternion_scalar = attitude_quaternion(1);
attitude_quaternion_vector = attitude_quaternion(2:4);
q = attitude_quaternion;
q0 = attitude_quaternion_scalar;
qbar = attitude_quaternion_vector(:);

%% Quaternion Kinematic Matrix
    % The attitude quaternion, q, has dynamics "qdot = E(q)*omega".  The
    % following "attitude_quaternion_kinematic_matrix" is the "E(q)" in 
    % that equation.
attitude_quaternion_kinematic_matrix = 1/2*[-qbar'; q0*eye(3) + ...
    hat_map(qbar)];

%% Rotation Matrix
    % The rotation matrix can be recalculated from the attitude quaternion.
    % Given the quaternion kinematics, "qdot = E(q)*omega", the rotation
    % matrix FROM THE BODY FRAME TO THE INERTIAL FRAME, R (i.e. v_inertial 
    % = R*v_body), is calculated
    %   R  = (q0^2 - qbar'*qbar)*eye(3) + 2*qbar*qbar' + 2*q0*hatmap(qbar)
    % and the rotation matrix FROM THE INERTIAL FRAM TO THE BODY FRAME, R'
    % (i.e. v_body = R'*v_inertial), is calculated
    %   R' = (q0^2 - qbar'*qbar)*eye(3) + 2*qbar*qbar' - 2*q0*hatmap(qbar)
    % There is a serious lack of clarity in many standard texts about which
    % formula should be used.  The above equations are correct.  This code
    % will implement the first of these two equations.  R' will be obtained
    % by transposing R.  The equation below is an equivalent method of
    % calculating R.  
attitude_quaternion_co_kinematic_matrix = 1/2*[-qbar'; q0*eye(3) - ...
    hat_map(qbar)];
attitude_matrix = 4*attitude_quaternion_co_kinematic_matrix'*...
    attitude_quaternion_kinematic_matrix;

%% Vestigal
% A previous version of this code integrated the attitude matrix.  The
% following line is an artifact of that version and may be useful in the
% future.  
    % attitude_quaternion = attitude_matrix_to_quaternion(attitude_matrix,previous_attitude_quaternion);
end
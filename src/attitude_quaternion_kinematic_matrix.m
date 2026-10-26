function attitude_quaternion_kinematic_matrix_value = attitude_quaternion_kinematic_matrix(attitude_quaternion)
%ATTITUDE_QUATERNION_KINEMATIC_MATRIX  The matrix E(q) in qdot = E(q)*omega_B.
%   E = ATTITUDE_QUATERNION_KINEMATIC_MATRIX(q) returns the 4x3 matrix that
%   maps body angular velocity to the time derivative of the body-to-inertial
%   unit quaternion,
%       qdot = E(q)*omega_B,
%   for the scalar-first convention q = [q0; qbar], with q0 the scalar part
%   and qbar the 3x1 vector part. In closed form,
%       E(q) = (1/2)*[      -qbar';
%                            q0*I + hat(qbar) ].
%
%   Convention: The rotation matrices in this code are body-to-inertial: R
%   (not R') maps vectors from the body frame to the inertial frame,
%   v_inertial = R*v_body. The sign of the q0*I + hat(qbar) block is tied to
%   that choice. Many texts present the E corresponding to the
%   inertial-to-body convention, in which the sign of the hat(qbar) term
%   differs; the form above is the one consistent with
%   ATTITUDE_QUATERNION_TO_MATRIX, which builds R from the same pair of
%   matrices.
%
%   The quaternion is assumed to be a unit quaternion. No normalization is
%   performed here: MULTIROTOR_FLOW_MAP normalizes q_IB before calling, and
%   MULTIROTOR_MANIFOLD_MAP renormalizes after every accepted step.
%
%   Inputs:
%       attitude_quaternion  4x1 unit quaternion, scalar first
%
%   Outputs:
%       attitude_quaternion_kinematic_matrix_value  4x3 matrix E(q)
%
%   Requires on the MATLAB path: hat_map.
%
%   Called by MULTIROTOR_FLOW_MAP.

%% Quaternion Components

%% Quaternion Components
    % This code uses the scalar first convention for qiaternions and
    % abbreviates the scalar part as q0 and the vector part as qbar.  
attitude_quaternion_scalar = attitude_quaternion(1);
attitude_quaternion_vector = attitude_quaternion(2:4);
q0 = attitude_quaternion_scalar;
qbar = attitude_quaternion_vector(:);

%% Quaternion Kinematic Matrix
    % The attitude quaternion, q, has dynamics "qdot = E(q)*omega".  The
    % following "attitude_quaternion_kinematic_matrix" is the "E(q)" in 
    % that equation.
attitude_quaternion_kinematic_matrix_value = 1/2*[-qbar'; q0*eye(3) + ...
    hat_map(qbar)];

end
function attitude_quaternion = attitude_matrix_to_quaternion(attitude_matrix,previous_attitude_quaternion)
% ATTITUDE_MATRIX_TO_QUATERNION Computes the attitude quaternion
%  corresponding to a given attitude matrix.  This algorithm uses the scalar
%  first convention.  
%
%  Note, since two (antipolar) unit quaternions correspond to each rotation 
%  matrix, the sign is ambiguous.  This algorithm computes the unit
%  quaternion with a positive scalar part and then checks which of the two
%  corresponding unit quaternions is closer to the previous unit
%  quaternion.  If no previous quaternion is supplied, the algorithm takes
%  the unit quaternion with the positive real component.  
%
%  Note, the common formulation of the formula converting from an attitude
%  matrix to an attitude quaternion assumes that the attitude matrix
%  rotates vectors from the INERTIAL FRAME to the BODY FRAME.  THIS CODE
%  USES THE OPPOSITE CONVENTION.  OUR BASIC ROTATION MATRIX, R (not R'),
%  ROTATES VECTORS FROM THE BODY FRAME TO THE INERTIAL FRAME.  Accordingly,
%  the following code has been designed to match our convention.  
%% Check for supplied previous quaternion...
if ~exist('previous_attitude_quaternion','var') || isempty(previous_attitude_quaternion)
    previous_attitude_quaternion = [1;0;0;0];
end

%% Compute scalar part...
unsigned_scalar = real(1/2*sqrt(1+trace(attitude_matrix)));

%% Compute vector part...
% The following computation was modified from the standard algorithm in the
% literature to account for the fact that the attitude matrix being taken
% as an input for this function is the one that rotates vectors from the
% body frame to the inertial frame (i.e. vi = attitude_matrix*vb).  The
% modification was to simply switch the sign of the vector part of the
% quaternion. 
unsigned_vector = 1/4/unsigned_scalar*...
    [   attitude_matrix(3,2)-attitude_matrix(2,3);
        attitude_matrix(1,3)-attitude_matrix(3,1);
        attitude_matrix(2,1)-attitude_matrix(1,2)   ];
if unsigned_scalar == 0
    [V,D] = eig(attitude_matrix); % The rotation axis is always an eigenvector of R
    d = diag(D);
    % The corresponding eigenvalue is always 1, but there may be numerical errors
    [~,flag] = min(abs(d-1));
    a = V(:,flag);
    unsigned_vector = a/norm(a);
end
%% Assemble unsigned quaternion...
unsigned_quaternion = [unsigned_scalar; unsigned_vector];

%% Determine sign of quaternion...
if unsigned_quaternion.'*previous_attitude_quaternion >= 0
    quaternion_sign = 1;
else
    quaternion_sign = -1;
end

%% Compute signed quaternion...
attitude_quaternion = quaternion_sign*unsigned_quaternion;
attitude_quaternion = attitude_quaternion/norm(attitude_quaternion);
end
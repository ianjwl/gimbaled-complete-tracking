function V = multirotor_lyapunov_function(tj,xctj,xdtj,parameters)
%% Initialization partition.
% Parse inputs...
t = tj(1);
j = tj(2);
vehicle_state = xctj;
controller_state = xdtj;

% Call user input function...
user_input = user_input_function(t,vehicle_state,parameters);

% Parse controller states...
IL_controller_state = controller_state(1:7,1);
reference_state = controller_state(12:29,1);
OL_controller_state = controller_state(30:end,1);

%% Reference model partition.
% Recall reference output...
C_ref = parameters.C_ref;
D_ref = parameters.D_ref;
reference_output = C_ref*reference_state + D_ref*user_input;

% Handle the case when reference is overridden...
if parameters.override_reference
    [x_user,dx_user,ddx_user,dddx_user,ddddx_user,~] = alternate_reference_function(t);
    reference_output(1:15,1) = [x_user;dx_user;ddx_user;dddx_user;ddddx_user];
end

%% Outer-loop controller partition.
% Load parameter values...
m = parameters.m;
g = parameters.g;
d_OL = [0;0;m*g];
A_OL = parameters.A_OL;
B_OL = parameters.B_OL;
C_OL = parameters.C_OL;
invCAB_OL = parameters.invCAB_OL;
K_P = parameters.K_P_OL;
K_D = parameters.K_D_OL;
Gamma = parameters.Gamma_OL;
C_jet_2 = [C_OL; C_OL*A_OL];
C_jet_3 = [C_OL; C_OL*A_OL; C_OL*A_OL^2];
KPD1 = [K_P K_D eye(3)];


% Separate vehicle states...
x_I = vehicle_state(1:3); % Inertial frame cartesian position
v_I = vehicle_state(4:6); % Velocity of the center of mass relative to the inertial frame resolved in the inertial frame

% Separate reference states...
y_jet_01 = reference_output(1:6);
y_jet_02 = reference_output(1:9);
y_jet_13 = reference_output(4:12);

% Parse controller states...
alpha_hat_OL = OL_controller_state;

% Load outer-loop control states...
x_OL = [x_I; v_I];

% Compute outer-loop control signal (target control)...
a_star = -d_OL + invCAB_OL*(alpha_hat_OL - KPD1*C_jet_3*x_OL + KPD1*y_jet_02);

% Prepare to compute ddt_x_OL...
ddt_x_OL_MATRIX = [A_OL*x_OL + B_OL*d_OL   B_OL];
    % ddt_x_OL = ddt_x_OL_MATRIX*[1; a];

% Prepare to compute ddt_a_star...
ddt_a_star_MATRIX = invCAB_OL*(...
    [Gamma*(y_jet_01 - C_jet_2*x_OL) + KPD1*y_jet_13 zeros(3,3)] ...
    -KPD1*C_jet_3*ddt_x_OL_MATRIX ...
    );
    % ddt_a_star = ddt_a_star_MATRIX*[1; a];

%% Inner-loop controller partition.

% Load parameter values...
r = parameters.r;
rho = parameters.rho;
base = (1-rho)/2;
alpha_2pi = base/4;
alpha_0 = max([0 (3*rho-1)/2]);
alpha_1 = min([2*rho (1+rho)/2]);
delta0 = parameters.delta0;
delta1 = parameters.delta1;
varrho = parameters.varrho;
Delta_a = parameters.Delta_a;
Delta_ddt_n_star = parameters.Delta_ddt_n_star;
k_a = parameters.k_a;
k_eta = parameters.k_eta;
a0 = parameters.a0;
a1 = parameters.a1;
k_mu = parameters.k_mu;
p_omega = parameters.p_omega;
R_CB = [unit_vector(1,3) -unit_vector(2,3) -unit_vector(3,3)];

% Separate vehicle states...
q_IB = vehicle_state(7:10); q_IB = q_IB/norm(q_IB); % Unit quaternion corresponding to R_IB
omega_B = vehicle_state(11:13); % Angular velocity of the body frame relative to the inertial frame resolved in the body frame
% Separate reference output...
nu = reference_output(16:18);
% Separate controller states...
mu_star = IL_controller_state(1:3,1);
qeta_lp = IL_controller_state(4:7,1);

% Preprocess attitude information...
R_IB = attitude_quaternion_to_matrix(q_IB);
R_IC = R_IB*R_CB';
R = R_IC;
lambda = R(:,3);
omega_C = R_CB*omega_B;
omega = omega_C;

% Calculate scalar thrust command...
a_ = norm(a_star);

% Select behavior based on whether the thrust command is below the threshold...
if a_ <= min([alpha_2pi*r a0])
    % In this case, V_a == 0 and omega_star == 0.  Hence the Lyapunov
    % function's value is
    V = p_omega/2*norm(omega)^2;
    return;
end

% Calculate target thrust direction...
lambda_star = a_star/a_;
% Calculate attitude error theta = acos(lambda_star'*lambda) with
% better numerical conditioning...
theta = atan2(norm(cross(lambda_star, lambda)), lambda_star'*lambda);
% Calculate theta_star...
alpha = a_/r;
theta_star = theta_star_function(alpha,rho);
% Calculate kappa_o...
kappa_o = rho + intssat(alpha,alpha_0,alpha_1);
% Calculate kappa...
ssat_kappa = ssat(theta-theta_star,delta0,delta1);
kappa = kappa_o*(1-ssat_kappa) + rho*ssat_kappa;
% Calculate u1...
u1 = r*kappa;

% Calculate MRP attitude error target...
if theta_star >= 2*pi
    n_star = inf; % <------ MATLAB returns "tan(pi/2) = 1.6331e+16" not "tan(pi/2) = inf"... Therefore, I'm setting the infinity manually!
else
    n_star = tan(theta_star/4);
end

% Set and correct mu_star...
Pi_lambda_star = (lambda_star*lambda_star');
mu_star = (eye(3) - Pi_lambda_star)*mu_star;
if norm(mu_star) == 0; error('Error: mu_star was parallel to lambda_star.'); end
mu_star = mu_star/norm(mu_star);
% Calculate Omega_mu...
Omega_mu = k_mu*cross(lambda_star,mu_star)'*nu;
% Calculate a...
a = lambda*u1;
% Calculate ddt_a_star...
ddt_a_star = ddt_a_star_MATRIX*[1; a];
% Calculate ddt_a_...
ddt_a_ = lambda_star'*ddt_a_star;
% Calculate ddt_alpha...
ddt_alpha = ddt_a_/r;
% Calculate R_star...
R_star = [mu_star cross(lambda_star,mu_star) lambda_star];
% Calculate R_tilde...
R_tilde = R_star'*R;
% Calculate qeta...
qeta = attitude_matrix_to_quaternion(R_tilde,qeta_lp); % <--- Always in the qeta_lp hemisphere
% Calculate eta...
eta = quat2mrp(qeta); % Default value which matches previous value...

% Calculate Omega_a...
Omega_a = k_a*intssat(a_*(norm(eta)-n_star),0,Delta_a)/ssat(a_*(norm(eta)-n_star),0,Delta_a);
% Calculate ddt_theta_star...
grad_theta_star = dtheta_star(alpha,rho);
ddt_theta_star = grad_theta_star*ddt_alpha;
% Calculate ddt_n_star...
ddt_n_star = sec(theta_star/4)^2*ddt_theta_star/4;
% Calculate omega_star_term_1...
omega_star_term_1 = -4*eta/norm(eta)/(1+norm(eta)^2)*Omega_a/a_*ssat(norm(eta)/n_star,varrho,1);
% Calculate omega_star_term_2...
omega_star_term_2 = R'*lambda_star*Omega_mu*ssat(norm(eta)/n_star,varrho,1);
% Calculate omega_star_term_3...
omega_star_term_3 = (R'*hat_map(lambda_star) - 4*eta/norm(eta)/(1+norm(eta)^2)*(norm(eta)-n_star)*lambda_star')*ddt_a_star/a_*ssat(norm(eta)/n_star,varrho,1);
% Calculate omega_star_term_4...
omega_star_term_4 = 4*eta/norm(eta)/(1+norm(eta)^2)*(ddt_n_star*(1-ssat(ddt_n_star,0,Delta_ddt_n_star)))*ssat(norm(eta)/n_star,varrho,1);
% Calculate omega_star_term_5...
omega_star_term_5 = -4*k_eta*eta/(1 + norm(eta)^2)*ssat(a_,a0,a1);

% Calculate omega_star...
if norm(eta)/n_star > varrho
    if isnan(omega_star_term_1)
        omega_star_term_1 = zeros(3,1); 
    end
    omega_star = ...
          omega_star_term_1 ...
        + omega_star_term_2 ...
        + omega_star_term_3 ...
        + omega_star_term_4 ...
        + omega_star_term_5;
else
    omega_star = omega_star_term_5;
end

% Calculate V_a...
V_a = intssat(a_*(norm(eta)-n_star),0,Delta_a);
% Calculate Lyapunov value...
V = V_a + p_omega/2*norm(omega - omega_star)^2;

end
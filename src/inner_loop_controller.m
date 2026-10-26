function [u_baseline,IL_controller_state_tpdt,outputs] = inner_loop_controller(t,dt,vehicle_state,reference_output,IL_controller_state,a_star,ddt_a_star_MATRIX,ddt2_a_star_MATRIX,parameters,outputs)
%INNER_LOOP_CONTROLLER  Gimbal control law: collective thrust and torque.
%   [u_baseline,xd_IL_tpdt,outputs] = INNER_LOOP_CONTROLLER(t,dt, ...
%       vehicle_state,reference_output,xd_IL,a_star,ddt_a_star_MATRIX, ...
%       ddt2_a_star_MATRIX,parameters,outputs)
%   synthesizes the collective thrust u1 and the body torque u2 that drive
%   the gimbaled actuator toward the ball A(t) = a*(t) + r*B^3, and returns
%   them stacked as u_baseline = [u1; u2]. u2 is a torque, not an angular
%   acceleration: the desired control-frame angular acceleration is rotated
%   into the body frame and passed through the rigid-body relation
%       u2 = In*alpha_B - omega_B x (In*omega_B).
%   u_baseline is therefore directly compatible with the vehicle dynamics in
%   MULTIROTOR_FLOW_MAP.
%
%   This is the backstepping design of Sections III-A-4 and III-A-6. The
%   function proceeds in two stages.
%
%   Stage 1: thrust. With alpha = ||a*||/r,
%       u1 = r*( kappa_circ(alpha)*(1 - zeta(theta - theta*,delta0,delta1))
%                + rho*zeta(theta - theta*,delta0,delta1) )
%   i.e. eq. (u1_prescription) with kappa_circ from eq. (defn_of_kappa_circ),
%   kappa from eq. (defn_of_kappa), and theta* from THETA_STAR_FUNCTION. The
%   radially asymmetric factor attenuates u1 when lambda and lambda* are
%   misaligned. When a* = 0 the law reduces to u1 = U1min, which is the
%   paper's extension to the free-fall singularity where theta is undefined.
%
%   Stage 2: angular acceleration. The MRP error eta is formed from the
%   error rotation R_tilde = R_star'*R, with the target frame
%   R_star = [mu_star, cross(lambda_star,mu_star), lambda_star] and mu_star
%   advanced by the auxiliary-direction objective eq.
%   (necessary_for_mu_str_dot) with Omega_mu = k_mu*cross(lambda_star,mu_star)'nu.
%   The target rate omega* is the five-term expression of eq.
%   (defn_of_omega_star): the Omega_a dissipation term, the Omega_mu and
%   ddt_a* feedforward terms, the ddt_n* containment term, and the k_eta
%   trim. Backstepping then gives u2 = ddt_omega_star - (1/2)k_omega
%   (omega - omega*) plus the eta-proportional cross term, i.e.
%   eq. (u2_prescription).
%
%   Regime threshold: alpha_2pi = (1-rho)/8 is the top of the 2*pi plateau
%   of theta*, so a_ < alpha_2pi*r implies n_star = inf and N_star = R^3. The
%   weaker condition a_ < a0 is used for the eta-independent branch as well,
%   because zeta(a_;a0,a1) is zero below a0, so the k_eta trim term also
%   vanishes there and the whole omega* collapses to zero. Testing
%   theta_star >= 2*pi instead would be wrong on [a0, alpha_2pi*r]; see the
%   in-file comment.
%
%   Containment gate: Whenever ||eta||/n_star <= varrho, only the k_eta trim
%   term of omega* is retained, since the other four are gated by
%   zeta(||eta||/n_star;varrho,1). This gate is what makes omega* continuous
%   through the free-fall singularity: it switches off the terms
%   proportional to ddt_a*/||a*|| before ||a*|| reaches zero.
%
%   Controller state: The suffix "_lp" abbreviates "last pass", i.e. the
%   value carried over from the previous call to this controller.
%       xd_IL(1:3)   mu_star, the auxiliary direction
%       xd_IL(4:7)   qeta_lp, the last-pass MRP-hemisphere memory
%                    quaternion, passed to ATTITUDE_MATRIX_TO_QUATERNION so
%                    that the extracted qeta stays in a consistent
%                    hemisphere from step to step
%   The returned state is advanced with a single forward-Euler increment on
%   mu_star, followed by renormalization, and with qeta_lp_tpdt = qeta.
%   qeta_lp is written unconditionally, so the hemisphere memory never goes
%   stale.
%
%   Note that the local mu_star is overwritten in the a* = 0 branch, where 
%   it is set directly to nu/||nu|| rather than projected onto the plane
%   orthogonal to lambda_star, since lambda_star is undefined there.
%
%   Naming caution: Inside this function, u1 is the collective thrust and u2
%   is the body torque. The name u2_alpha is used for the intermediate
%   desired angular acceleration, before conversion to a torque.
%
%   Inputs:
%       t, dt              current time and flow step length; t is unused
%       vehicle_state      layout as in MULTIROTOR_FLOW_MAP; uses
%                          vehicle_state(7:13), the attitude and body rate
%       reference_output   21x1 output of REFERENCE_MODEL; uses 16:21
%       xd_IL              7x1 inner-loop state, layout above
%       a_star             3x1 target actuation from OUTER_LOOP_CONTROLLER
%       ddt_MATRIX         3x4 coefficient matrix for ddt_a*
%       ddt2_MATRIX        3x7 coefficient matrix for ddt2_a*
%       parameters         struct of system and controller parameters; uses
%                          In, U1min, r, rho, delta0, delta1, varrho,
%                          Delta_a, Delta_ddt_n_star, k_a, k_eta, a0, a1,
%                          k_mu, k_omega, p_omega
%       outputs            diagnostics struct, threaded through
%
%   Outputs:
%       u_baseline         4x1 [u1; u2], collective thrust then body torque
%       xd_IL_tpdt         7x1 inner-loop state at t + dt
%       outputs            diagnostics struct, extended with roughly forty
%                          entries, including theta, theta_star, the five
%                          omega* terms and their time derivatives, Omega_a,
%                          V, ddt_V, ddt_V_a, eta, R_tilde, the a*
%                          derivative chain, ddt_omega_R_star (the time
%                          derivative of the target-frame angular velocity),
%                          and more. Every entry that is only defined in one
%                          branch is guarded with exist(...,"var"), so the
%                          struct has the same fields on every call; entries
%                          not computed on a given call are NaN-filled.
%
%   Requires on the MATLAB path:
%       attitude_quaternion_to_matrix, attitude_matrix_to_quaternion,
%       quat2mrp, mrp_kinematic_matrix, hat_map, unit_vector,
%       theta_star_function, intssat, ssat, dssat, dtheta_star,
%       ddtheta_star.
%
%   Called by CONTROL_LAW.
%
%   IMPORTANT: MULTIROTOR_JUMP_GUARD and MULTIROTOR_LYAPUNOV_FUNCTION both
%   rebuild this omega* construction and the formula for V inline, rather
%   than calling this function, so that they can evaluate the guard and the
%   Lyapunov value without running the whole cascade. Any change to the
%   control law here must be mirrored in both files, or the switching test
%   and the step-size control will silently disagree with the control that
%   is actually applied. 

% Load parameter values...
In = parameters.In;
U1min = parameters.U1min;
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
k_omega = parameters.k_omega;
p_omega = parameters.p_omega;
R_CB = [unit_vector(1,3) -unit_vector(2,3) -unit_vector(3,3)];

% Separate vehicle states...
q_IB = vehicle_state(7:10); q_IB = q_IB/norm(q_IB); % Unit quaternion corresponding to R_IB
omega_B = vehicle_state(11:13); % Angular velocity of the body frame relative to the inertial frame resolved in the body frame

% Separate reference output...
nu = reference_output(16:18);
ddt_nu = reference_output(19:21);

% Separate controller states...
mu_star = IL_controller_state(1:3,1); 
qeta_lp = IL_controller_state(4:7,1); 


%% Inner-loop Controller: Propulsion Calculation

% Preprocess attitude information...
R_IB = attitude_quaternion_to_matrix(q_IB);
R_IC = R_IB*R_CB';
R = R_IC;
lambda = R(:,3);
% mu = R(:,1);
omega_C = R_CB*omega_B;
omega = omega_C;

% Calculate scalar thrust command...
a_ = norm(a_star); 

% Select behavior based on whether the thrust command is zero...
if a_ == 0
    
    % Set attitude error target...
    theta_star = 2*pi;

    % Set value of u1...
    u1 = U1min;

else % if a_ ~= 0...

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

end



%% Inner-loop Controller: Angular Acceleration Calculation

% Calculate MRP attitude error target...
if theta_star >= 2*pi
    n_star = inf; % <------ MATLAB returns "tan(pi/2) = 1.6331e+16" not "tan(pi/2) = inf"... Therefore, I'm setting the infinity manually!
else
    n_star = tan(theta_star/4);
end

% Select behavior based on whether the thrust command is zero...
if a_ == 0

    % Set desired (control frame) angular acceleration...
    u2_alpha = -k_omega/2*omega;

    % Set error quaternion value for persistency...
    qeta = qeta_lp; % The value of eta is actually undefined when a_ == 0, but this choice ensures qeta_lp is well defined on the next evaluation

    % Set value and rate of pointing vector target for persistency...
    if norm(nu) > 0
        mu_star = nu/norm(nu);
        ddt_mu_star = (eye(3) - (nu*nu')/norm(nu)^2)*ddt_nu/norm(nu);
    else
        ddt_mu_star = zeros(3,1);
    end

    % Set Lyapunov function value for persistency...
    V = 1/2*p_omega*norm(omega)^2;
    % V_a = 0;

    % Set Lyapunov rate...
    ddt_V = p_omega*omega'*u2_alpha;
    % ddt_V_a = 0;

else % if a_ ~= 0...

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

    % Calculate ddt_lambda_star...
    ddt_lambda_star = (eye(3) - Pi_lambda_star)*ddt_a_star/a_;

    % Calculate ddt_mu_star...
    ddt_mu_star = Omega_mu*cross(lambda_star,mu_star) - lambda_star*(mu_star'*ddt_lambda_star);

    % Calculate R_star...
    R_star = [mu_star cross(lambda_star,mu_star) lambda_star];

    % Calculate R_tilde...
    R_tilde = R_star'*R;

    % Calculate qeta...
    qeta = attitude_matrix_to_quaternion(R_tilde,qeta_lp); % <--- Always in the qeta_lp hemisphere

    % Calculate eta...
    eta = quat2mrp(qeta); % Default value which matches previous value...

    % Select behavior based on whether the control law is independent of eta...
    if a_ < min([alpha_2pi*r a0]) % Then N_star = \RR^3 and the control law is independent of eta!      <------ Used to be "if isinf(n_star)" but MATLAB returns "tan(pi/2) = 1.6331e+16" not "tan(pi/2) = inf"... Also used to be "if theta_star >= 2*pi", but this is wrong if a_ is not less than a0!

        % Set omega_star...
        omega_star = zeros(3,1);

        % Set desired (control frame) angular acceleration...
        u2_alpha = -k_omega*omega/2;

        % Set Lyapunov function value for persistency...
        V = 1/2*p_omega*norm(omega)^2;
        % V_a = 0;

        % Set Lyapunov rate...
        ddt_V = p_omega*omega'*u2_alpha;
        % ddt_V_a = 0;

    else % if a_ >= min([alpha_2pi*r a0]), then the control law does depend on eta!

        % Calculate Omega_a...
        Omega_a = k_a*intssat(a_*(norm(eta)-n_star),0,Delta_a)/ssat(a_*(norm(eta)-n_star),0,Delta_a);
        Omega_a_is_valid = true;
        if isnan(Omega_a) || isinf(Omega_a)
            Omega_a = 0;
            Omega_a_is_valid = false;
        end
        
        % Calculate omega_R_star...
        omega_R_star = [
            -cross(lambda_star,mu_star)'*ddt_lambda_star
            mu_star'*ddt_lambda_star
            Omega_mu
        ];

        % Calculate grad_theta_star and grad2_theta_star...
        grad_theta_star = dtheta_star(alpha,rho);
        grad2_theta_star = ddtheta_star(alpha,rho);

        % Calculate ddt_theta_star...
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
        
        % Calculate Lyapunov value V...
        V = V_a + p_omega/2*norm(omega - omega_star)^2;

        % Calculate ddt_lambda...
        ddt_lambda = cross( R*omega, lambda );

        % Calculate ddt_theta = ddt_theta = -(lambda_star'*ddt_lambda + lambda'*ddt_lambda_star)/sqrt(1-(lambda'*lambda_star)^2) with a better numerical approach...
        if abs(lambda'*lambda_star) ~= 1
            ddt_theta = (cross(lambda,ddt_lambda) - cross(lambda_star,ddt_lambda_star))'*(cross(lambda_star,lambda)/norm(cross(lambda_star,lambda)));
        else
            ddt_theta = 0; % In this case, ddt_theta should be unused... (True if theta_star is designed correctly...)
        end

        % Calculate grad_kappa_o...
        grad_kappa_o = ssat(alpha,alpha_0,alpha_1);

        % Calculate ddarg1_kappa...
        gradarg1_kappa = grad_kappa_o*(1-ssat(theta-theta_star,delta0,delta1)) ...
            + (kappa_o - rho)*dssat(theta-theta_star,delta0,delta1)*grad_theta_star;

        % Calculate ddarg2_kappa...
        gradarg2_kappa = -(kappa_o - rho)*dssat(theta-theta_star,delta0,delta1);
        
        % Calculate ddt_kappa...
        if gradarg1_kappa ~= 0
            temp1 = gradarg1_kappa*ddt_alpha;
        else
            temp1 = 0;
        end
        if gradarg2_kappa ~= 0
            temp2 = gradarg2_kappa*ddt_theta;
        else
            temp2 = 0;
        end
        ddt_kappa = temp1 + temp2;

        % Calculate ddt_a...
        ddt_a = r*ddt_kappa*lambda + ...
                r*kappa*ddt_lambda;

        % Calculate ddt2_a_star...
        ddt2_a_star = ddt2_a_star_MATRIX*[1; a; ddt_a];

        % Calculate ddt2_a_...
        ddt2_a_ = lambda_star'*ddt2_a_star + ddt_lambda_star'*ddt_a_star;

        % Calculate ddt2_alpha...
        ddt2_alpha = ddt2_a_/r;

        % Calculate ddt2_lambda_star...
        ddt2_lambda_star = 1/a_*(ddt2_a_star - 2*ddt_lambda_star*ddt_a_ - lambda_star*ddt2_a_);

        % Calculate ddt_Omega_mu...
        ddt_Omega_mu = k_mu*cross(ddt_lambda_star,mu_star)'*nu + ...
            k_mu*cross(lambda_star,ddt_mu_star)'*nu + ...
            k_mu*cross(lambda_star,mu_star)'*ddt_nu;
        
        % Calculate ddt_omega_R_star... 
        ddt_omega_R_star = [
            -cross(lambda_star,mu_star)'*ddt2_lambda_star - cross(lambda_star,ddt_mu_star)'*ddt_lambda_star
            mu_star'*ddt2_lambda_star + ddt_mu_star'*ddt_lambda_star
            ddt_Omega_mu
        ];

        % Calculate omega_tilde...
        omega_tilde = omega - R_tilde'*omega_R_star;

        % Calculate ddt_R_tilde...
        ddt_R_tilde = R_tilde*hat_map(omega_tilde);

        % Calculate ddt_eta...
        ddt_eta = mrp_kinematic_matrix(eta)*omega_tilde;

        % Calculate ddt_Omega_a...
        if Omega_a_is_valid
            ddt_Omega_a = (k_a - Omega_a*...
                dssat(a_*(norm(eta)-n_star),0,Delta_a)/ssat(a_*(norm(eta)-n_star),0,Delta_a))*...
                    (ddt_a_*(norm(eta)-n_star) + a_*((eta/norm(eta))'*ddt_eta - ddt_n_star));
        else
            ddt_Omega_a = 0;
        end

        % Calculate ddt2_theta_star...
        ddt2_theta_star = grad2_theta_star*ddt_alpha^2 + grad_theta_star*ddt2_alpha;

        % Calculate ddt2_n_star...
        ddt2_n_star = 1/8*sec(theta_star/4)^2*(2*ddt2_theta_star + ddt_theta_star^2*tan(theta_star/4));
        % ddt2_n_star = ddt_n_star*(ddt2_theta_star/ddt_theta_star + n_star/2);
        
        % Calculate eta_ratio_term and ddt_eta_ratio_term...
        eta_ratio_term = 4*eta/norm(eta)/(1+norm(eta)^2);
        ddt_eta_ratio_term = (-8*(eta*eta')/norm(eta)/(1+norm(eta)^2)^2 ...
            + 4/(1+norm(eta)^2)*(eye(3) ...
            - (eta*eta')/norm(eta)^2)/norm(eta))*ddt_eta;

        % Calculate saturation_ratio_term and ddt_saturation_ratio_term...
        saturation_ratio_term = ssat(norm(eta)/n_star,varrho,1);
        ddt_saturation_ratio_term = dssat(norm(eta)/n_star,varrho,1)*(...
            eta'*ddt_eta/norm(eta)/n_star-norm(eta)/n_star^2*ddt_n_star);

        % Calculate ddt_omega_star_term_1...
        ddt_omega_star_term_1 = ...
            - ddt_eta_ratio_term*Omega_a/a_*saturation_ratio_term ...
            - eta_ratio_term*ddt_Omega_a/a_*saturation_ratio_term ...
            + eta_ratio_term*Omega_a/a_^2*ddt_a_*saturation_ratio_term ...
            - eta_ratio_term*Omega_a/a_*ddt_saturation_ratio_term;

        % Calculate ddt_omega_star_term_2...
        ddt_omega_star_term_2 = ...
              (R*hat_map(omega))'*lambda_star*Omega_mu*saturation_ratio_term ...
            + R'*ddt_lambda_star*Omega_mu*saturation_ratio_term ...
            + R'*lambda_star*ddt_Omega_mu*saturation_ratio_term ...
            + R'*lambda_star*Omega_mu*ddt_saturation_ratio_term;

        % Calculate ddt_omega_star_term_3...
        ddt_omega_star_term_3 = ...
              (R*hat_map(omega))'*hat_map(lambda_star)*ddt_a_star/a_*saturation_ratio_term ...
            + R'*hat_map(ddt_lambda_star)*ddt_a_star/a_*saturation_ratio_term ...
            + R'*hat_map(lambda_star)*ddt2_a_star/a_*saturation_ratio_term ...
            - R'*hat_map(lambda_star)*ddt_a_star/a_^2*ddt_a_*saturation_ratio_term ...
            + R'*hat_map(lambda_star)*ddt_a_star/a_*ddt_saturation_ratio_term ...
            - ddt_eta_ratio_term*(norm(eta)-n_star)*lambda_star'*ddt_a_star/a_*saturation_ratio_term ...
            - eta_ratio_term*(eta'*ddt_eta/norm(eta)-ddt_n_star)*lambda_star'*ddt_a_star/a_*saturation_ratio_term ...
            - eta_ratio_term*(norm(eta)-n_star)*ddt_lambda_star'*ddt_a_star/a_*saturation_ratio_term ...
            - eta_ratio_term*(norm(eta)-n_star)*lambda_star'*ddt2_a_star/a_*saturation_ratio_term ...
            + eta_ratio_term*(norm(eta)-n_star)*lambda_star'*ddt_a_star/a_^2*ddt_a_*saturation_ratio_term ...
            - eta_ratio_term*(norm(eta)-n_star)*lambda_star'*ddt_a_star/a_*ddt_saturation_ratio_term;

        % Calculate ddt_omega_star_term_4...
        ddt_omega_star_term_4 = ...
              ddt_eta_ratio_term*ddt_n_star*(1-ssat(ddt_n_star,0,Delta_ddt_n_star))*saturation_ratio_term ...
            + eta_ratio_term*ddt2_n_star*(1-ssat(ddt_n_star,0,Delta_ddt_n_star))*saturation_ratio_term ...
            + eta_ratio_term*ddt_n_star*(-dssat(ddt_n_star,0,Delta_ddt_n_star)*ddt2_n_star)*saturation_ratio_term ...
            + eta_ratio_term*ddt_n_star*(1-ssat(ddt_n_star,0,Delta_ddt_n_star))*ddt_saturation_ratio_term;

        % Calculate ddt_omega_star_term_5...
        ddt_omega_star_term_5 = ...
            -4*k_eta*ddt_eta/(1+norm(eta)^2)*ssat(a_,a0,a1) ...
            +8*k_eta*(eta*eta')*ddt_eta/(1+norm(eta)^2)^2*ssat(a_,a0,a1) ...
            -4*k_eta*eta/(1+norm(eta)^2)*dssat(a_,a0,a1)*ddt_a_;

        % Calculate ddt_omega_star... 
        if norm(eta)/n_star > varrho
            ddt_omega_star = ...
                  ddt_omega_star_term_1 ...
                + ddt_omega_star_term_2 ...
                + ddt_omega_star_term_3 ...
                + ddt_omega_star_term_4 ...
                + ddt_omega_star_term_5;
        else
            ddt_omega_star = ddt_omega_star_term_5;
        end

        % Calculate the angular acceleration cross-term...
        cross_term = - 1/p_omega*(1+norm(eta)^2)/4*(eta/norm(eta))*a_*...
                    ssat(a_*(norm(eta)-n_star),0,Delta_a);

        % Calculate desired (control frame) angular acceleration... 
        if norm(eta) > n_star
            u2_alpha = ddt_omega_star - k_omega/2*(omega - omega_star) + cross_term;
        else
            u2_alpha = ddt_omega_star - k_omega/2*(omega - omega_star);       
        end

        % Calculate theoretical value of ddt_V...
        ddt_V = -(  Omega_a + ddt_n_star*a_*ssat(ddt_n_star,0,Delta_ddt_n_star) + k_eta*norm(eta)*a_*ssat(a_,a0,a1)  )*ssat(a_*(norm(eta)-n_star),0,Delta_a) ...
            - k_omega*p_omega/2*norm(omega - omega_star)^2;
        ddt_V_a = -(  Omega_a + ddt_n_star*a_*ssat(ddt_n_star,0,Delta_ddt_n_star) + k_eta*norm(eta)*a_*ssat(a_,a0,a1)  )*ssat(a_*(norm(eta)-n_star),0,Delta_a) ...
            + (1+norm(eta)^2)/4*a_*ssat(a_*(norm(eta)-n_star),0,Delta_a)/norm(eta)*eta'*(omega-omega_star);

    end


end

% Convert desired control frame angular acceleration to body frame...
alpha_C = u2_alpha;
alpha_B = R_CB'*alpha_C;

% Convert desired (body frame) angular acceleration to torque... 
u2 = In*alpha_B - cross(In*omega_B,omega_B);

% Assemble inner-loop control output...
u_baseline = [u1; u2];

% Compute mu_star_tpdt...
mu_star_tpdt = mu_star + ddt_mu_star*dt;
mu_star_tpdt = mu_star_tpdt/norm(mu_star_tpdt);

% Compute qeta_lp_tpdt...
qeta_lp_tpdt = qeta;

% Assemble inner-loop control state increment...
IL_controller_state_tpdt = [mu_star_tpdt; qeta_lp_tpdt];


%% Assemble and Organize Outputs
outputs.u_baseline = u_baseline;
if ~exist("eta","var"); eta = nan(3,1); end
outputs.eta = eta;
outputs.mu_star = mu_star;
outputs.ddt_mu_star = ddt_mu_star;
outputs.mu_star_tpdt = mu_star_tpdt;
outputs.lyapunov_value = V;
if ~exist("theta","var"); theta = nan(1,1); end
outputs.theta = theta;
if ~exist("ddt_theta","var"); ddt_theta = nan(1,1); end
outputs.ddt_theta = ddt_theta;
outputs.theta_star = theta_star;
if ~exist("ddt_theta_star","var"); ddt_theta_star = 0; end
outputs.ddt_theta_star = ddt_theta_star;
if ~exist("saturation_ratio_term","var"); saturation_ratio_term = nan; end
outputs.saturation_ratio_term = saturation_ratio_term;
if ~exist("ddt_saturation_ratio_term","var"); ddt_saturation_ratio_term = nan; end
outputs.ddt_saturation_ratio_term = ddt_saturation_ratio_term;
if ~exist("ddt2_theta_star","var"); ddt2_theta_star = 0; end
outputs.ddt2_theta_star = ddt2_theta_star;
if ~exist("omega_star","var"); omega_star = zeros(3,1); end
outputs.omega_star = omega_star;
outputs.omega_star_B = R_CB'*omega_star;
if ~exist("switch_boolean","var"); switch_boolean = false; end
outputs.switching_variable = switch_boolean;
if ~exist("ddt_omega_star","var"); ddt_omega_star = zeros(3,1); end
outputs.ddt_omega_star = ddt_omega_star;
outputs.n_star = n_star;
if ~exist("ddt_n_star","var"); ddt_n_star = 0; end
outputs.ddt_n_star = ddt_n_star;
if ~exist("ddt2_n_star","var"); ddt2_n_star = 0; end
outputs.ddt2_n_star = ddt2_n_star;
if ~exist("Omega_a","var"); Omega_a = 0; end
outputs.Omega_a = Omega_a;
if ~exist("ddt_Omega_a","var"); ddt_Omega_a = 0; end
outputs.ddt_Omega_a = ddt_Omega_a;
if ~exist("ddt_eta","var"); ddt_eta = nan(3,1); end
outputs.ddt_eta = ddt_eta;
if ~exist("R_tilde","var"); R_tilde = zeros(3,3); end
outputs.R_tilde_col = R_tilde(:);
if ~exist("ddt_R_tilde","var"); ddt_R_tilde = zeros(3,3); end
outputs.ddt_R_tilde_col = ddt_R_tilde(:);
if ~exist("a_","var"); a_ = 0; end
outputs.a_ = a_;
if ~exist("ddt_a_","var"); ddt_a_ = 0; end
outputs.ddt_a_ = ddt_a_;
if ~exist("ddt2_a_","var"); ddt2_a_ = 0; end
outputs.ddt2_a_ = ddt2_a_;
outputs.a_star = a_star;
if ~exist("ddt_a_star","var"); ddt_a_star = ddt_a_star_MATRIX*[1; lambda*u1]; end
outputs.ddt_a_star = ddt_a_star;
if ~exist("ddt2_a_star","var"); ddt2_a_star = nan(3,1); end
outputs.ddt2_a_star = ddt2_a_star;
if ~exist("lambda","var"); lambda = nan(3,1); end
outputs.lambda = lambda;
if ~exist("ddt_lambda","var"); ddt_lambda = nan(3,1); end
outputs.ddt_lambda = ddt_lambda;
if ~exist("lambda_star","var"); lambda_star = nan(3,1); end
outputs.lambda_star = lambda_star;
if ~exist("ddt_lambda_star","var"); ddt_lambda_star = nan(3,1); end
outputs.ddt_lambda_star = ddt_lambda_star;
if ~exist("ddt_V","var"); ddt_V = nan(1,1); end
outputs.ddt_V = ddt_V;
outputs.ddt_omega_B = (In\(u2 + cross(In*omega_B,omega_B)));
if ~exist("omega_star_term_1","var")
    omega_star_term_1 = nan(3,1); 
    omega_star_term_2 = nan(3,1); 
    omega_star_term_3 = nan(3,1); 
    omega_star_term_4 = nan(3,1); 
    omega_star_term_5 = nan(3,1); 
end
outputs.omega_star_term_1 = omega_star_term_1;
outputs.omega_star_term_2 = omega_star_term_2;
outputs.omega_star_term_3 = omega_star_term_3;
outputs.omega_star_term_4 = omega_star_term_4;
outputs.omega_star_term_5 = omega_star_term_5;
if ~exist("ddt_omega_star_term_1","var")
    ddt_omega_star_term_1 = nan(3,1); 
    ddt_omega_star_term_2 = nan(3,1); 
    ddt_omega_star_term_3 = nan(3,1); 
    ddt_omega_star_term_4 = nan(3,1); 
    ddt_omega_star_term_5 = nan(3,1); 
end
outputs.ddt_omega_star_term_1 = ddt_omega_star_term_1;
outputs.ddt_omega_star_term_2 = ddt_omega_star_term_2;
outputs.ddt_omega_star_term_3 = ddt_omega_star_term_3;
outputs.ddt_omega_star_term_4 = ddt_omega_star_term_4;
outputs.ddt_omega_star_term_5 = ddt_omega_star_term_5;
if ~exist("cross_term","var"); cross_term = nan(3,1); end
outputs.cross_term = cross_term;
if ~exist("ddt_omega_R_star","var"); ddt_omega_R_star = nan(3,1); end
outputs.ddt_omega_R_star = ddt_omega_R_star;
if ~exist("ddt_V_a","var"); ddt_V_a = nan(1,1); end
outputs.ddt_V_a = ddt_V_a;


end
function [xctjp1,xdtjp1] = multirotor_jump_map(tj,xctj,xdtj,parameters)

% Most states do not jump...
xctjp1 = xctj;
xdtjp1 = xdtj;

% The only state that jumps is the quaternion corresponding to eta_lp...
qeta_lp_tj = xdtj(4:7,1);

% The value it jumps to is the equivalent quaternion on the other hemisphere...
qeta_lp_tjp1 = -qeta_lp_tj;

% Update outputs...
xdtjp1(4:7,1) = qeta_lp_tjp1;

end
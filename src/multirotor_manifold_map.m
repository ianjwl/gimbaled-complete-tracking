function [xctj_out,xdtj_out] = multirotor_manifold_map(tj,xctj_in,xdtj_in,parameters)

% Ensure unit vectors are unit vectors...
xctj_out = [xctj_in(1:6,1); xctj_in(7:10)/norm(xctj_in(7:10)); xctj_in(11:13)];
xdtj_out = [xdtj_in(1:3,1)/norm(xdtj_in(1:3,1)); xdtj_in(4:7,1)/norm(xdtj_in(4:7,1)); xdtj_in(8:end,1)];

end
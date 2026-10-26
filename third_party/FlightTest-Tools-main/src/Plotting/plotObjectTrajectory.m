function plotObjectTrajectory(translations,rotations,scaleFactor,numObjects,stl,indices)
%plotObjectTrajectory
%
% Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
%
% This function plots the trajectory of a rigid body. It displays the path
% as well as a number of 3D models depicting attitude.
%
% Requirements: UAV Toolbox
%
% Inputs:
%
%   translations  See documentation for plotTransforms.m
%
%   rotations     See documentation for plotTransforms.m
%
%   scaleFactor   A scale factor that determines how large the 3D object
%                 defined by the input 'stl' is displayed.
%
%   numObjects    The number of 3D objects to plot along the path.
%
%   stl           The filename of an STL file that defines the 3D object.
%
%

% input arguments error checking
if (nargin<5)||(nargin>6)
    error('Incorrect number of inputs.');
end
if (length(translations)~=length(rotations))
    error('Number of samples is not consistent.');
end

% color of path
Ectoplasm = [0 255/255 165/255]; % [213/255 94/255 0]; % [0 144/255 178/255]; % [230/255 159/255 0]; % [0 158/255 115/255]; % 
pathColor = Ectoplasm;

% Check which type of rotation is supplied
if isa(rotations, "so3")
    % ROTATIONS is an SO object
    % Extract quaternion rotations
    rotQCol = quaternion(rotm(rotations), "rotmat", "point");
elseif isa(rotations, "so2")
    % Extract 2D rotation and change to 3D
    rotQCol = so2ToQuatCol(rotations);
else
    % ROTATIONS is a quaternion
    % validate rotation using internal quaternion validation utility
    robotics.internal.validation.validateQuaternion(rotations, 'plotTransforms', 'rotations');

    if (isa(rotations, 'quaternion') && size(rotations, 2) ~= 1)
        error(message('shared_robotics:robotcore:plotTransforms:RotationMustBeColumn'));
    end

    % convert orientations to quaternions if necessary
    if ~isa(rotations, 'quaternion')
        rotQCol = quaternion(rotations);
    else
        rotQCol = rotations;
    end
end

% NED --> ENU
ENU = [translations(:,2), translations(:,1), -translations(:,3)];
x = ENU(:,1);
y = ENU(:,2);
z = ENU(:,3);
R_BS = so3([pi/2 0 pi],"eul");
R_IB = so3(rotQCol);
R_IS = R_IB*R_BS;
R_EI = so3([pi/2 0 pi],"eul");
R_ES = R_EI*R_IS;
rotations = R_ES;

% number of samples
N = length(translations);

% plot the object path
hold on
plot3(x,y,z,'linewidth',1.5,'Color',pathColor)

% plot the objects
if nargin == 5
    samplingFactor = floor(N/(numObjects-1));
    indices = 1:samplingFactor:N;
end
for k = indices
    plotTransforms(ENU(k,:),rotations(k,:),'MeshFilePath',stl,'FrameSize',scaleFactor,'MeshColor',[0.75 0.75 0.75])
end
hold off
grid on

% Axis labels and legend
xlabel('East [m]','FontSize',10)
ylabel('North [m]','FontSize',10)
zlabel('Altitude [m]','FontSize',10)

% set the axes limits
limits = [min(x),max(x),min(y),max(y),min(z),max(z)];
limits = limits + scaleFactor*[-1,1,-1,1,-1,1];
axis(limits);

% set lighting, view, aspect ratio, etc.
grid on
view(35,25);
daspect([1 1 1]);
lightangle(35,60)

end

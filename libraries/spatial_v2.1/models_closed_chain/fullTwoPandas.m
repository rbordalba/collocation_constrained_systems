function  robot = fullTwoPandas

persistent memory;

if ~isempty(memory)
  robot = memory;
  return
end

%initial x0

%PARAMETERS
robot.NB = 14; %excluding base
robot.parent = [0 1 2 3 4 5 6 0 8 9 10 11 12 13];
robot.jtype = {'Rz','Rz','Rz','Rz','Rz','Rz','Rz','Rz','Rz','Rz','Rz','Rz','Rz','Rz'};
robot.loop{1} = [-1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1];
r = 0.159;
xoffset = 0.32;
d = 0.25; % vertical offset to the center of end effector
robot.Xloop{1} =  (rotz(-pi/2)*xlt([2*r,0,0])*rotz(pi/2));
% end effector information
robot.ee = [1 1 1 1 1 1 1 0 0 0 0 0 0 0];
robot.Xee =  plux(rz(0),[r,0,d])*plux(rz(pi/2),[0,0,0]);

%transformation from frame i-1 to i
% IT IS THE INVERSE OF THE CUIK TRANSFORM.. ROTX IS ROTX^-1 (SIGN% CHANGES)
%PANDA1
robot.Xtree{1} = plux(eye(3),[0,0,0.333]);
m(1) = 3.2512;
J{1} = [0.00043075 -6.9984e-07 -1.3145e-06;-6.9984e-07   0.0010312    0.010811;-1.3145e-06    0.010811   0.0046025];
CM{1} = [1.14e-05   -0.031266   -0.069373];

robot.Xtree{2} =  plux(rx(-pi/2),[0,0,0]);
m(2) = 3.2832;
J{2} = [0.00025946  8.2912e-07 -2.1162e-06;8.2912e-07   0.0046712   -0.011019;-2.1162e-06   -0.011019  0.00084776];
CM{2} = [-1.32e-05   -0.070322    0.031178];

robot.Xtree{3} =  plux(rx(pi/2),[0,-0.316,0]);  
m(3) = 2.5439;
J{3} = [0.0028132   0.0012348   -0.006994;0.0012348   0.0010618  -0.0038729;-0.006994  -0.0038722   0.0016002];
CM{3} = [0.044338    0.024919   -0.038137];

robot.Xtree{4} = plux(rx(pi/2),[0.0825,0,0]); 
m(4) = 2.5942;
J{4} = [0.0028013  -0.0011542  -0.0040781;-0.0011542   0.0029125   0.0040474;-0.0040781   0.0040474   0.0022384];
CM{4} = [-0.03855    0.039526    0.024716];

robot.Xtree{5} =  plux(rx(-pi/2),[-0.0825,0.384,0]); 
m(5) = 3.7365;
J{5} = [-0.018092 -1.9284e-05 -1.7342e-05;-1.9284e-05   -0.015541   -0.022195;-1.7342e-05   -0.022195   0.0029348];
CM{5} = [-6.37e-05    0.038412    -0.10997];

robot.Xtree{6} = plux(rx(pi/2),[0,0,0]); 
m(6) = 1.5682;
J{6} = [0.0023987  0.00037631  0.00081502;0.00037631 -0.00026964  9.8923e-05;0.00081502  9.8923e-05   0.0005599];
CM{6} = [0.051002   0.0069327    0.006169];

robot.Xtree{7} = plux(eye(3),[0,0,0.107])*plux(rx(pi/2),[0.088,0,0]);  
m(7) = 0.48744;
J{7} = [-0.0026303 -9.2185e-05  0.00038372;-9.2185e-05  -0.0026307   0.0003836;0.00038372   0.0003836  0.00063893];
CM{7} = [0.010361     0.01036    0.079108];     

%PANDA2
robot.Xtree{8} = plux(eye(3),[0,0,0.333])*plux(rz(0*pi),[xoffset,0,0]);  
m(8) = 3.2512;
J{8} = [0.00043075 -6.9984e-07 -1.3145e-06;-6.9984e-07   0.0010312    0.010811;-1.3145e-06    0.010811   0.0046025];
CM{8} = [1.14e-05   -0.031266   -0.069373];

robot.Xtree{9} = plux(rx(-pi/2),[0,0,0]);
m(9) = 3.2832;
J{9} = [0.00025946  8.2912e-07 -2.1162e-06;8.2912e-07   0.0046712   -0.011019;-2.1162e-06   -0.011019  0.00084776];
CM{9} = [-1.32e-05   -0.070322    0.031178];

robot.Xtree{10} =  plux(rx(pi/2),[0,-0.316,0]);  
m(10) = 2.5439;
J{10} = [0.0028132   0.0012348   -0.006994;0.0012348   0.0010618  -0.0038729;-0.006994  -0.0038722   0.0016002];
CM{10} = [0.044338    0.024919   -0.038137];

robot.Xtree{11} = plux(rx(pi/2),[0.0825,0,0]); 
m(11) = 2.5942;
J{11} = [0.0028013  -0.0011542  -0.0040781;-0.0011542   0.0029125   0.0040474;-0.0040781   0.0040474   0.0022384];
CM{11} = [-0.03855    0.039526    0.024716];

robot.Xtree{12} =  plux(rx(-pi/2),[-0.0825,0.384,0]); 
m(12) = 3.7365;
J{12} = [-0.018092 -1.9284e-05 -1.7342e-05;-1.9284e-05   -0.015541   -0.022195;-1.7342e-05   -0.022195   0.0029348];
CM{12} = [-6.37e-05    0.038412    -0.10997];

robot.Xtree{13} = plux(rx(pi/2),[0,0,0]); 
m(13) = 1.5682;
J{13} = [0.0023987  0.00037631  0.00081502;0.00037631 -0.00026964  9.8923e-05;0.00081502  9.8923e-05   0.0005599];
CM{13} = [0.051002   0.0069327    0.006169];

robot.Xtree{14} = plux(eye(3),[0,0,0.107])*plux(rx(pi/2),[0.088,0,0]);  
m(14) = 0.48744;
J{14} = [-0.0026303 -9.2185e-05  0.00038372;-9.2185e-05  -0.0026307   0.0003836;0.00038372   0.0003836  0.00063893];
CM{14} = [0.010361     0.01036    0.079108];    

for i=1:robot.NB
    %Inertia  (m, CoM, Inertia_z)
    robot.I{i}=mcI(m(i),CM{i},J{i});
end

%%%%robot.constraint = @constraint; % constraint-imposing function
%%%%robot.CartesianConstraint = @CartesianConstraint; % End Effector constraint

%gravity vector
robot.gravity = [0; 0; -9.81]; 

p1 = 0.07/5;
p2 = 0.04/5;
R = 0.1/5;
robot.appearance.base = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{1} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{2} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{3} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{4} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{5} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{6} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{7} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{8} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{9} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{10} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{11} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{12} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{13} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{14} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], R };

robot.n_u = 14;
robot.n_e = 12; %includes velocity constraints
% Projection into action space  tau = B*u;
robot.B = eye(14);

robot.friction = ones(robot.NB,1)*0.2;
% Limits on action and states
robot.umin = [-87,-87,-87,-87,-12,-12,-12,-87,-87,-87,-87,-12,-12,-12]';
robot.umax = [87,87,87,87,12,12,12,87,87,87,87,12,12,12]';
%robot.xmin = [-pi/2-0.2; -2; -2*pi; -pi/2; -2*pi; -2*pi; -2*pi;      -pi/2-0.2; -2;  -2*pi; -pi/2; -2*pi; -2*pi; -2*pi; repmat(-25,robot.NB,1)];
robot.xmin = [-pi/2-0.1; -2; -2*pi; -2*pi; -2*pi; -2*pi; -pi-0.1;      -pi/2-0.1; -2;  -2*pi; -2*pi; -2*pi; -2*pi; -pi-0.1; repmat(-25,robot.NB,1)];
%robot.xmax = [-pi/2+0.2; +2; 2*pi; pi/2; 2*pi; -pi/2; 2*pi;           -pi/2+0.2; +2; 2*pi; pi/2; 2*pi; -pi/2; 2*pi; repmat(25,robot.NB,1)];
robot.xmax = [-pi/2+0.1; +2; 2*pi; 2*pi; 2*pi; -pi/2; -pi+0.1;           -pi/2+0.1; +2; 2*pi; 2*pi; 2*pi; -pi/2; -pi+0.1; repmat(25,robot.NB,1)];
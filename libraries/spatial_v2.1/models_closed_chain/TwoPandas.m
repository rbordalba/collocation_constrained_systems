function  robot = TwoPandas

persistent memory;

if ~isempty(memory)
  robot = memory;
  return
end

%initial x0
robot.init = [-1.9169257028484532 0.0083206323833330 -1.2699529738031643 0.0850060524085256 2.4864700296939706 -1.8967530714947185 0.0007657812527465 -1.1673973819015515 0.0875276352800123 2.4037799208376196 0 0 0 0 0 0 0 0 0 0]';
robot.init = [-0.381166216164 0.034795155495 -0.562323792509 0.045483288179 -3.032679421770 -0.344227656391 0.027426010656 -0.432067219814 0.051417573098 -3.126124141637 -0.447291544242 -0.090498712273 -1.050102881807 0.191045014530 -0.924673816870 -0.591121922729 -0.127291546734 -1.244143025447 0.222274761915 -0.875139827217];
robot.init = [-0.381166216164 0.034795155495 -0.562323792509 0.045483288179 -3.032679421770 -0.344227656391 0.027426010656 -0.432067219814 0.051417573098 -3.126124141637 0 0 0 0 0 0 0 0 0 0 ]';

q0 = robot.init;

%PARAMETERS
robot.NB = 10; %excluding base
robot.parent = [0 1 2 3 4 0 6 7 8 9];
robot.jtype = {'Rx','Rx','Rx','Rx','Rx','Rx','Rx','Rx','Rx','Rx'};
robot.loop{1} = [-1 -1 -1 -1 -1 1 1 1 1 1]; % loop information (-1 parent branch, 1 successor branch)

% CONSTANTS
r = 0.159;
d = 0.25;
x_offset = 0.32;
mass_load = 12;
radius_load = r;
height_load=0.4;

% JOINT TRANSFORMS
%transformation from frame i-1 to i
%PANDA1
% from link0 to link1 (fixed), from link1 to link 2 (Featherstone inverse
% transform)
robot.Xtree{1} = rotx(pi)*rotz(pi/2)*plux(rz(-pi/2),[0,0,0.333]);
robot.T{1} = makehgtform('translate',[0 0 0.333],'zrotate',-pi/2)*...
    makehgtform('zrotate',pi/2)*makehgtform('xrotate',pi);
m(1) = 3.2832;
J{1} = [0.00025946  8.2912e-07 -2.1162e-06;8.2912e-07   0.0046712   -0.011019;-2.1162e-06   -0.011019  0.00084776];
CM{1} = [-1.32e-05   -0.070322    0.031178];
X2L{1} = rotx(pi/2)*rotz(pi/2); % X from joint to link

robot.Xtree{2} = plux(rz(-pi/2),[0,-0.316,0])*rotx(pi/2)*rotz(pi/2);
robot.T{2} = makehgtform('zrotate',pi/2)*makehgtform('xrotate',pi/2)*...
             makehgtform('translate',[0 -0.316 0],'zrotate',-pi/2); 
m(2) = 2.5439;
J{2} = [0.0028132   0.0012348   -0.006994;0.0012348   0.0010618  -0.0038729;-0.006994  -0.0038722   0.0016002];
CM{2} = [0.044338    0.024919   -0.038137];
X2L{2} = rotx(pi/2)*rotz(pi/2); % X from joint to link

robot.Xtree{3} = plux(rz(-pi/2),[0.0825,0,0])*rotx(pi/2)*rotz(pi/2);
robot.T{3} = makehgtform('zrotate',pi/2)*makehgtform('xrotate',pi/2)*...
makehgtform('translate',[0.0825 0 0],'zrotate',-pi/2);
m(3) = 2.5942;
J{3} = [0.0028013  -0.0011542  -0.0040781;-0.0011542   0.0029125   0.0040474;-0.0040781   0.0040474   0.0022384];
CM{3} = [-0.03855    0.039526    0.024716];
X2L{3} = rotz(pi/2)*roty(pi/2); % X from joint to link

robot.Xtree{4} = rotz(pi/2)*plux(ry(pi),[-0.0825,0.384,0])*rotz(pi/2)*roty(pi/2);
robot.T{4} =  makehgtform('yrotate',pi/2)*makehgtform('zrotate',pi/2)*...
    makehgtform('translate',[-0.0825 0.384 0],'yrotate',pi)*makehgtform('zrotate',pi/2);
m(4) = 3.7365;
J{4} = [-0.018092 -1.9284e-05 -1.7342e-05;-1.9284e-05   -0.015541   -0.022195;-1.7342e-05   -0.022195   0.0029348];
CM{4} = [-6.37e-05    0.038412    -0.10997];
X2L{4} = rotz(pi/2)*roty(pi/2); % X from joint to link

robot.Xtree{5} = rotz(-pi/2)*rotz(pi/2)*roty(pi/2);
robot.T{5} = makehgtform('yrotate',pi/2)*makehgtform('zrotate',pi/2)*...
makehgtform('zrotate',-pi/2); 
m(5) = 1.5682;
J{5} = [0.0023987  0.00037631  0.00081502;0.00037631 -0.00026964  9.8923e-05;0.00081502  9.8923e-05   0.0005599];
CM{5} = [0.051002   0.0069327    0.006169];
X2L{5} = rotx(pi/2)*rotz(pi/2);%rotx(pi/2)*rotz(pi/2); % X from joint to link

% fixed joint 
% transformation from succeesor to predecesor in the 0-DOF joint (p^X_s)
robot.Tloop{1} = makehgtform('zrotate',pi/2)*makehgtform('xrotate',pi/2)*...
makehgtform('translate',[0.088 0 0],'xrotate',pi/2)*makehgtform('zrotate',-pi)*...
    makehgtform('translate',[0 0 0.107]);
robot.Tloop{1} = robot.Tloop{1}*makehgtform('zrotate',pi/2)*makehgtform('translate',[r 0 d]);
robot.Tloop{2} = makehgtform('zrotate',pi/2)*makehgtform('xrotate',pi/2)*...
    makehgtform('translate',[0.088 0 0],'xrotate',pi/2)*makehgtform('zrotate',-pi)*...
    makehgtform('translate',[0 0 0.107]);
robot.Tloop{2} = robot.Tloop{2}*makehgtform('zrotate',pi/2)*makehgtform('translate',[-r 0 d]);

% from 6 to hand (inverse transform than 
robot.Xloop{1} = xlt([0,0,0.107])*rotz(-pi)*rotx(pi/2)*xlt([0.088,0,0])*rotx(pi/2)*rotz(pi/2);
% form hand to object
robot.Xloop{1} = xlt([r,0,d])*rotz(pi/2)*robot.Xloop{1};
% from object to 6'
robot.Xloop{1} = inv(xlt([-r,0,d])*rotz(pi/2)*xlt([0,0,0.107])*rotz(-pi)*rotx(pi/2)*xlt([0.088,0,0])*rotx(pi/2)*rotz(pi/2))*robot.Xloop{1};

%PANDA2
robot.Xtree{6} = rotx(pi)*rotz(pi/2)*plux(rz(-pi/2),[0,0,0.333])*xlt([x_offset,0,0]);
robot.T{6} = makehgtform('translate',[x_offset 0 0])*makehgtform('translate',[0 0 0.333],'zrotate',-pi/2)*...
    makehgtform('zrotate',pi/2)*makehgtform('xrotate',pi);
m(6) = 3.2832;
J{6} = [0.00025946  8.2912e-07 -2.1162e-06;8.2912e-07   0.0046712   -0.011019;-2.1162e-06   -0.011019  0.00084776];
CM{6} = [-1.32e-05   -0.070322    0.031178];
X2L{6} = rotx(pi/2)*rotz(pi/2); % X from joint to link

robot.Xtree{7} = plux(rz(-pi/2),[0,-0.316,0])*rotx(pi/2)*rotz(pi/2);
robot.T{7} = makehgtform('zrotate',pi/2)*makehgtform('xrotate',pi/2)*...
             makehgtform('translate',[0 -0.316 0],'zrotate',-pi/2); % correct!?
m(7) = 2.5439;
J{7} = [0.0028132   0.0012348   -0.006994;0.0012348   0.0010618  -0.0038729;-0.006994  -0.0038722   0.0016002];
CM{7} = [0.044338    0.024919   -0.038137];
X2L{7} = rotx(pi/2)*rotz(pi/2); % X from joint to link

robot.Xtree{8} = plux(rz(-pi/2),[0.0825,0,0])*rotx(pi/2)*rotz(pi/2);
robot.T{8} = makehgtform('zrotate',pi/2)*makehgtform('xrotate',pi/2)*...
makehgtform('translate',[0.0825 0 0],'zrotate',-pi/2);
m(8) = 2.5942;
J{8} = [0.0028013  -0.0011542  -0.0040781;-0.0011542   0.0029125   0.0040474;-0.0040781   0.0040474   0.0022384];
CM{8} = [-0.03855    0.039526    0.024716];
X2L{8} = rotz(pi/2)*roty(pi/2); % X from joint to link

robot.Xtree{9} = rotz(pi/2)*plux(ry(pi),[-0.0825,0.384,0])*rotz(pi/2)*roty(pi/2);
robot.T{9} =  makehgtform('yrotate',pi/2)*makehgtform('zrotate',pi/2)*...
    makehgtform('translate',[-0.0825 0.384 0],'yrotate',pi)*makehgtform('zrotate',pi/2);
m(9) = 3.7365;
J{9} = [-0.018092 -1.9284e-05 -1.7342e-05;-1.9284e-05   -0.015541   -0.022195;-1.7342e-05   -0.022195   0.0029348];
CM{9} = [-6.37e-05    0.038412    -0.10997];
X2L{9} = rotz(pi/2)*roty(pi/2); % X from joint to link

robot.Xtree{10} = rotz(-pi/2)*rotz(pi/2)*roty(pi/2);
robot.T{10} = makehgtform('yrotate',pi/2)*makehgtform('zrotate',pi/2)*...
makehgtform('zrotate',-pi/2); 
m(10) = 1.5682;
J{10} = [0.0023987  0.00037631  0.00081502;0.00037631 -0.00026964  9.8923e-05;0.00081502  9.8923e-05   0.0005599];
CM{10} = [0.051002   0.0069327    0.006169];
X2L{10} = rotx(pi/2)*rotz(pi/2);%rotx(pi/2)*rotz(pi/2); % X from joint to link

%extended inertia for other fixed bodies (Featherstone Handbook)
Xlink7 = rotz(-pi)*rotx(pi/2)*xlt([0.088,0,0])*rotx(pi/2)*rotz(pi/2);
Xlink8 = xlt([0,0,0.107])*Xlink7;
Xhand = Xlink8;
Xobject1 = xlt([r,0,d])*rotz(pi/2)*Xhand;
Xobject2 = xlt([-r,0,d])*rotz(pi/2)*Xhand;
m7 = 0.48744;
CM7 = [0.010361,0.01036,0.079108]';
J7 = [-0.0026303,-9.2185e-05,0.00038372;-9.2185e-05,-0.0026307,0.0003836;0.00038372,0.0003836,0.00063893];
I7 = Xlink7'*mcI(m7,CM7,J7)*Xlink7;
mh = 0.53538;
CMh = [-1.67e-5,0.0013237,0.027468]';
Jh = [0.0016372,4.2807e-08,-3.3302e-07; 4.2807e-08,-0.00011003,1.8822e-05;-3.3302e-07,1.8822e-05,0.0018777];
Ih = Xhand'*mcI(mh,CMh,Jh)*Xhand;
mobj = mass_load;
CMobj = [0,0,0]';
Jobj = [1/12*mass_load*(3*radius_load^2+height_load^2),0,0;
        0,1/12*mass_load*(3*radius_load^2+height_load^2),0;
        0,0,1/2*mass_load*radius_load^2];
Iobj = mcI(mobj,CMobj,Jobj)/2;
Iobj1 = Xobject1'*Iobj*Xobject1;
Iobj2 = Xobject2'*Iobj*Xobject2;

for i=1:robot.NB
    %Inertia  (m, CoM, Inertia_z)
    robot.I{i}=X2L{i}'*mcI(m(i),CM{i},J{i})*X2L{i};
    %robot.I{i}%
end

% add inertia of fixed bodies to bodies 5 and 10
robot.I{5} = robot.I{5}+I7+Ih+Iobj1;
robot.I{10} = robot.I{10}+I7+Ih+Iobj2;

%gravity vector
robot.gravity = [0; 0; 9.81];

% end effector info
robot.ee = [1 1 1 1 1 0 0 0 0 0]; % end effector route (at 3rd link coordinate system)
robot.cartesian = [4 5 6]'; % independent end-effector coordinates: x position, y position, z position
% transformation to end effector
robot.Xee = plux(eye(3),[0.0,0,-0.4]); % transform to end effector from ee link coord system


p1 = 0.07/5;
p2 = 0.04/5;
R = 0.1/5;
robot.appearance.base = ...
{ 'cyl', [0 0 -p1; 0 0 p1], 1.5*R };
robot.appearance.body{1} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], 1.5*R };
robot.appearance.body{2} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], 1.5*R };
robot.appearance.body{3} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], 1.5*R };
robot.appearance.body{4} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], 1.5*R };
robot.appearance.body{5} = ...
{ 'cyl', [0 0 -p1; 0 0 p1], 1.5*R };
robot.appearance.body{6} = ...
{ 'cyl', [0 0 -p2; 0 0 p2], R };
robot.appearance.body{7} = ...
{ 'cyl', [0 0 -p2; 0 0 p2], R };
robot.appearance.body{8} = ...
{ 'cyl', [0 0 -p2; 0 0 p2], R };
robot.appearance.body{9} = ...
{ 'cyl', [0 0 -p2; 0 0 p2], R };
robot.appearance.body{10} = ...
{ 'cyl', [0 0 -p2; 0 0 p2], R };

robot.n_u = 10;
robot.n_e = 6*2;
% Projection into action space  tau = B*u;
robot.B = eye(10);

robot.friction = ones(robot.NB,1)*0.2;

% Limits on action and states
robot.umin = [-87,-87,-87,-12,-12,-87,-87,-87,-12,-12]';
robot.umax = [87,87,87,12,12,87,87,87,12,12]';
robot.xmin = [-2; -2*pi; -2*pi; -2*pi; -2*pi; -2; -2*pi; -2*pi; -2*pi; -2*pi; repmat(-25,robot.NB,1)];
robot.xmax = [+2; +2*pi; +2*pi; +2*pi; +2*pi; +2; +2*pi; +2*pi; +2*pi; +2*pi;  repmat(25,robot.NB,1)];


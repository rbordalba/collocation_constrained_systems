function  robot = Scherbot2(p)

% FiveBars robot built by Beta Robots using Featherstone library and
% Ricards mods.

persistent memory;

if ~isempty(memory)
  robot = memory;
  return
end

%%%%%%%%%%%%%%%%%%%%%%% JOINTS DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
robot.NB = 5; %number of bodies/joints
robot.parent = [0 1 2 3 0];
robot.jtype = {'Rz','Rz','Rz','Rz','Rz'};
robot.loop{1} = [-1 -1 -1 -1 1]; % loop information (-1 parent branch, 1 successor branch)
% define axis that define the robot [rotation; cartesian] 
robot.planar = [3 4 5]'; % z rotation, x position, y position
robot.ee = [1 1 0 0 0]; % end effector route (at 3rd link coordinate system)
robot.cartesian = [4 5]'; % independent end-effector coordinates: x position, y position

%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
robot.gravity = [0; -9.81; 0]; % gravity vector

% End-effector
m_load = 0.5; % mass of the load 
L_load = 0.06; % radius of load cylinder 
% Lengths
L0 = 0.12; %base link length (measured in lab)
Lproximal = 0.2; % length of proximal link (measured in lab)
Ldistal = 0.15; % length of distal link (measured in lab)
L = [Lproximal Ldistal Ldistal Lproximal Lproximal]; %last link is reapeated with half its inertia
robot.L = [L0 L];
% center of mass length (not necessarely link length)
if nargin == 1
    Lcm = [p.Lcm(1); p.Lcm(2); L(3)-p.Lcm(2); L(4)-p.Lcm(1); p.Lcm(1)];
else
    Lcm = [0.11; 0.0746; Ldistal-0.0746; Lproximal-0.11; 0.11];
end

% Links masses
mdistal = 0.8925; %mass of distal links (estimated)
mproximal = 1.2093; % mass of proximal (estimated)
if nargin == 1
    m = [p.m(1) p.m(2) p.m(2) p.m(1)/2 p.m(1)/2];
else
    m = [mproximal mdistal mdistal mproximal/2 mproximal/2];
end
% link inertias
if nargin == 1
    Iz = [p.Iz(1); p.Iz(2); p.Iz(2); p.Iz(1)/2; p.Iz(1)/2];
else 
    Iz = [7 2.3 2.3 7/2 7/2]*1e-3;
end
%Center of mass, inertia 
for i=1:robot.NB
    CM{i} = [Lcm(i),0,0];
    J{i}=  diag([0 0 Iz(i)]);
    robot.I{i}=mcI(m(i),CM{i},J{i});
end

if m_load ~= 0
    CMl = [L(2),0,0]; % load is attached at the end of distal link
    Jl=  diag([2/5*m_load*L_load^2 2/5*m_load*L_load^2 2/5*m_load*L_load^2]); % inertia of sphere?
    % Combine inertia parameters of Load with link 2
    mt = m(2)+m_load;
    CMt = (CM{2}*m(2)+CMl*m_load)/mt;
    rl = norm(CMl-CMt);
    r2 = norm(CM{2}-CMt);
    Jt = J{2}+diag([0 0 m(2)*r2^2])+Jl+diag([0 0 m_load*rl^2]);
    robot.I{2}=mcI(mt,CMt,Jt);
end

% Actuation and friciton
if nargin == 1
    sc1 = p.sc(1); % estimated scaled factor Motor 1
    sc2 = p.sc(2); % estimated scaled factor Motor 2
    b1 = p.b(1);
    b2 = p.b(2);
    c1 = p.c(1);
    c2 = p.c(2);
else
    sc1 = 1; % estimated scaled factor Motor 1
    sc2 = 1; % estimated scaled factor Motor 2
    b1 = 0.0176;%1.76e-2;
    b2 = 0.0176; %6.9e-3;
    c1 = 0;%0.0665;
    c2 = 0;%0.0758; 
end
robot.n_u = 2; % number of actuators
robot.B = [sc1 0;0 0;0 0;0 0;0 sc2]; % Projection into action space  tau = B*u;
robot.n_e = 6; % number of position + velocity constraints
robot.friction = [b1 0 0 0 b2]'; % joint viscous friction
robot.Coulomb = 0*[c1 0 0 0 c2]'; % for Optimization 0

% Limits on action and states
robot.umin = repmat(-2,robot.n_u,1); % minimum allowed torque
robot.umax = repmat(2,robot.n_u,1); % maximum allowed torque
robot.xmin = [repmat(-10*pi,robot.NB,1); repmat(-30,robot.NB,1)];
robot.xmax = [repmat(10*pi,robot.NB,1);  repmat(30,robot.NB,1)];

%%%%%%%%%%%%%%%%%%%% GEOMETRIC TRANSFORMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%transformation from frame i-1 to i
robot.Xtree = {plux(rz(pi),[0,0,0]),plux(eye(3),[L(1),0,0]), plux(eye(3),[L(2),0,0]), plux(eye(3),[L(3),0,0]),plux(ry(pi)*rz(pi),[L0,0,0])};
% transformation from succeesor to predecesor in the 0-DOF joint (p^X_s)
robot.Xloop{1} =  plux(ry(-pi),[L(4),0,0]);
% transformation to end effector
robot.Xee = plux(eye(3),[L(2),0,0]); % transform to end effector from ee link coord system


%%%%%%%%%%%%%%%%%%% HARD-CODED CONSTRAINTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
robot.CartesianConstraint = @(x)  CartesianConstraint(robot,x(1:robot.NB),x(robot.NB+1:end)); % End Effector constraint
robot.FindState = @(posX,posY,velX,velY,mode) FindState(posX, posY, velX, velY, robot, mode);

%%%%%%%%%%%%%%% LINKS GRAPHICAL DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1 = 0.007;
p2 = 0.004;
R = 0.01;
robot.appearance.base = ...
{ 'box', [0 -p1 -p2; L0 p1 p2], ...
  'cyl', [0 0 -p1; 0 0 p1], 1.5*R };

robot.appearance.body{1} = ...
{ 'box', [0 -p1 -p2; L(1) p1 p2], ...
  'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{2} = ...
{ 'box', [0 -p1 -p2; L(2) p1 p2], ...
  'sphere', [L(2) 0 0], 2*R, ...
  'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{3} = ...
{ 'box', [0 -p1 -p2; L(3) p1 p2], ...
  'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{4} = ...
{ 'box', [0 -p1 -p2; L(4) p1 p2], ...
  'cyl', [0 0 -p1; 0 0 p1], R };
robot.appearance.body{5} = ...
{ 'box', [0 -p1 -p2; L(5) p1 p2], ...
  'cyl', [0 0 -p1; 0 0 p1], R };
             
function x0 = FindState(posX, posY, velX, velY, model, mode)

L0 = model.L(1);

n_q = model.NB;
Bx_l = 0; By_l = 0; % abs. coord of left base joint
Bx_r = L0; By_r = 0; % abs. coord of right base joint

% Find configuration given [x,y] coordinates
%mode = '-+';
%posX = 0.1; posY = 0.08354; %desired [x,y] position
[theta1,theta2] = InvKinLeg(Bx_l,By_l,posX,posY,model.L(2),model.L(3),mode(1));
[theta5,theta4] = InvKinLeg(Bx_r,By_r,posX,posY,model.L(2),model.L(3),mode(2));
theta3 = wrapToPi(theta2-theta4);
if isempty(theta3)
    error('The [x,y] position given is out of the workspace');
end

q = [theta1-pi,theta2-theta1,-theta3+pi,theta5-(theta2-theta3)]';
q = wrapToPi([q;wrapTo2Pi(2*pi-sum(q))]);
               
J1 = [model.L(2)*sin(q(1))+model.L(3)*sin(q(1)+q(2)) model.L(3)*sin(q(1)+q(2));
      -model.L(2)*cos(q(1))-model.L(3)*cos(q(1)+q(2)) -model.L(3)*cos(q(1)+q(2))];

J2 = [-model.L(3)*sin(q(5)+q(4))  -model.L(2)*sin(q(5))-model.L(3)*sin(q(4)+q(5));
      -model.L(3)*cos(q(5)+q(4)) -model.L(2)*cos(q(5))-model.L(3)*cos(q(4)+q(5))];

qDot = [J1 zeros(size(J1,1),1) zeros(size(J2)); ...
        -ones(1,length(q)); ...
        zeros(size(J1)) zeros(size(J1,1),1) J2]\[velX;velY;0;velX;velY];
x0 = [q;qDot];



function [theta1,theta2] = InvKinLeg(Bx,By,Px,Py,l1,l2,wm)
% Compute the inverse kinematics of a leg
%    Input:
%        Bx, By  Absolute coords. of base joint (Float)
%        Px, Py  Absolute coords. of the platform joint (Float)
%        l1, l2  Lengths of the proximal and distal links (Float)
%        wm      Working mode of the leg (Char, '+' or '-')
%
%    Output:
%        [theta1,theta2] The link angles relative to the absolute OX axis
%                        (Note that they are not wrapped to [-pi,pi])

% Cosine of the relative joint angle at the middle joint
dsquared = (Px-Bx)^2 + (Py-By)^2;
calpha2 = (dsquared - l1^2 - l2^2) / (2*l1*l2);

% Relative joint angle alpha2 at the middle joint
if wm == '-'
    alpha2 =  acos(calpha2);
else
    alpha2 = - acos(calpha2);
end


% Relative joint angle alpha1 at the proximal joint
A = l1+l2*cos(alpha2);
B = l2*sin(alpha2);
M = [A -B; B A];
result = M \ [Px-Bx;Py-By];
if ~isreal(result)
    theta1 = [];
    theta2 = [];
    return;
end  
alpha1 = atan2(result(2,1),result(1,1));

% Absolute joint angles (no need to wrapToPi() in our context)
theta1 = alpha1;
theta2 = alpha1 + alpha2;

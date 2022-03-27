function  [H,C,Phi,Phi_q,xi,Phi_q_dot] = HandC( model, q, qd, f_ext )

% HandC  Calculate coefficients of equation of motion
% [H,C]=HandC(model,q,qd,f_ext)  calculates the coefficients of the
% joint-space equation of motion, tau=H(q)qdd+C(d,qd,f_ext), where q, qd
% and qdd are the joint position, velocity and acceleration vectors, H is
% the joint-space inertia matrix, C is the vector of gravity,
% external-force and velocity-product terms, and tau is the joint force
% vector.  Algorithm: recursive Newton-Euler for C, and
% Composite-Rigid-Body for H.  f_ext is an optional argument specifying the
% external forces acting on the bodies.  It can be omitted if there are no
% external forces.  The format of f_ext is explained in the source code of
% apply_external_forces.
% Algorithms mainly from HandBook of Robotics 2nd Ed, Chapter Dynamics from
% Roy Featherstone. Table 2.9 for Composite Rigid body, Table 2.6 for
% Recursive Newton-Euler and 2.11 for computing loop-closure constraints.

% Forward pass of the recursive Newton-Euler (inverse dynamics with qDDot = 0)
% (Table 2.6 - Handbook of Robotics 2nd Ed.)
% Also forward pass to compute loop-closure constraints (Table 2.11 - Handbook of Robotics 2nd Ed.)
a_grav = get_gravity(model);
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    avp{i} = Xup{i} * -a_grav;
    
    % Table 2.11: Compute i^X_0
    Xa{i} = Xup{i}; 
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    avp{i} = Xup{i}*avp{model.parent(i)} + crm(v{i})*vJ;
    
    %NEW: Compute i^X_0
    Xa{i} = Xup{i} * Xa{model.parent(i)};
  end
  fvp{i} = model.I{i}*avp{i} + crf(v{i})*model.I{i}*v{i};
  
  % NEW: transform Si to 0 coordinates    0^X_i * S{i}
  J{i} = Xa{i}\S{i};
  Jdot{i} = Xa{i}\(crm(v{i})*S{i});
end

if nargin == 4
  fvp = apply_external_forces( model.parent, Xup, fvp, f_ext );
end

% Backward pass (Table 2.6 - Handbook of Robotics 2nd Ed.) (inverse dynamics with qDDot = 0)
% Backward pass composite inertia calculation (Table 2.9 - Handbook of Robotics 2nd Ed.)
IC = model.I;				
for i = model.NB:-1:1
  C{i} = S{i}' * fvp{i};
  if model.parent(i) ~= 0
    fvp{model.parent(i)} = fvp{model.parent(i)} + Xup{i}'*fvp{i};
    IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}'*IC{i}*Xup{i};
  end
end
C = vertcat(C{:});

if ~isnumeric(q)
    import casadi.*
    H = SX.zeros(model.NB,model.NB); % SX or MX depending on prev. defined
end			

% Forward pass composite inertia calculation (Table 2.9 - Handbook of Robotics 2nd Ed.)
for i = 1:model.NB
  fh = IC{i} * S{i};
  H(i,i) = S{i}' * fh;
  j = i;
  while model.parent(j) > 0
    fh = Xup{j}' * fh;
    j = model.parent(j);
    H(i,j) = S{j}' * fh;
    H(j,i) = H(i,j);
  end
end

J = horzcat(J{:});
Jdot = horzcat(Jdot{:});

%NEW
Phi = [];
Phi_q = [];
Phi_q_dot = [];
xi = [];

% Compute loop-closure constraints (Table 2.11 - Handbook of Robotics 2nd Ed.)
for i = 1:length(model.loop)
   p = max(find(model.loop{i}==-1));
   s = max(find(model.loop{i}==1));
%    % loop Transform
%    if isempty(s)
%        Xl = Xa{p}^-1*model.Xloop{i}^-1; % Featherstone in 0 coordinates 
%    else
%        Xl = Xa{p}^-1*model.Xloop{i}^-1*Xa{s}; % Featherstone in 0 coordinates
%    end
   mloop = repmat(model.loop{i},6,1);
   Phi_q = [Phi_q; mloop.*J]; % Featherstone in 0 coordinates
   Phi_q_dot = [Phi_q_dot; mloop.*Jdot];
  
   if isempty(s)
       xi = [xi;Xa{p}\avp{p}+a_grav];
       Phi= [Phi; X2vec(Xa{p}^-1*model.Xloop{i}^-1)]; % Featherstone in 0 coordinates
   else
       xi = [xi;Xa{p}\avp{p} - Xa{s}\avp{s}];
       Phi= [Phi; X2vec(Xa{p}^-1*model.Xloop{i}^-1*Xa{s})]; % Featherstone in 0 coordinates
   end
   
end

%Xa{p};
%Xa{s};

% Remove rows of zeros for planar mechanisms (only works for single-loop)
if isfield(model,'planar')
    dep = model.planar;
    Phi = Phi(dep);
    Phi_q = Phi_q(dep,:);
    Phi_q_dot = Phi_q_dot(dep,:);
    xi = xi(dep);
% else
%     dep = any(abs(Phi_q)>1e-10,2);
%     Phi = Phi(dep);
%     Phi_q = Phi_q(dep,:);
%     Phi_q_dot = Phi_q_dot(dep,:);
%     xi = xi(dep);
end

function vec = X2vec(X)

   R = X(1:3,1:3); rl = X(4:6,1:3);
   Sr = rl*R^-1; %obtain screw(r) in global coordinates

   %when finding config it works better
   vec = [X(2,3);
          X(3,1);
          X(1,2);
          Sr(2,3);
          Sr(3,1);
          Sr(1,2)];
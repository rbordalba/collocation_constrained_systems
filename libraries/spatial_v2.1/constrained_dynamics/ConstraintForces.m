function  [f_rel,f_abs] = ConstraintForces( model, x, u, f_ext)
% CONSTRAINTFORCES Computes constraint forces in closed-chain mechanisms
% ConstraintForces( model, x, u) compute the constraint forces transmited
% at each joint from child to parent in the inertia coordinate system
% (f_abs) and the link coordinate system (f_rel). 
% 
% The closed-chain system is modelled as a tree like mechanism, and a 
% generalized constraint force is computed by means of Lagrange multipliers
% that closes the closed-chain system. These Lagrange multipliers are the
% wrench applied at the cut link and is modelled as an external force to
% the tree like mechanism.
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu

n_q = model.NB;
q = x(1:n_q);
qd = x(n_q+1:end);

% Compute acceleration and lambdas
if nargin == 4
    [qdd, ~, lambda] = ForwardDynamics([],x,u,model,3,f_ext);
else
    [qdd, ~, lambda] = ForwardDynamics([],x,u,model,3);
end

% Compute external wrench closing the loop from lambdas s.t. Wrench = [Tau_loop; Force_loop];
if ~isnumeric(q)
    import casadi.*
    for i = 1:length(model.loop)  
        f_loop{i} = SX.zeros(6,1); 
    end
else
    for i = 1:length(model.loop)  
        f_loop{i} = zeros(6,1);
    end
end

if isfield(model,'planar')
    for i = 1:length(model.loop)  
        f_loop{i}(model.planar) = lambda(1:3)+(i-1)*3;
    end
else
    for i = 1:length(model.loop)  
        f_loop{i} = lambda((1:6)+(i-1)*6); 
    end
end
    
a_grav = get_gravity(model);

% Step 1 and 2 of Recursive Newton-Euler to compute acceleration of links
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  % velocity across joint i
  vJ = S{i}*qd(i);
  % Transformation from i-1 to i
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    % Modelling gravity as a virtual acc. at the base: Sec 5.3
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
  else
    % Velocity propagation: Eq.(5.7)-(5.14)-(5.18)
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    % acceleration propagation: Eq.(5.8)-(5.15)-(5.19)
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
  end
  % Compute net force acting on body i Eq.(5.9)
  f{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
end

% Step 3: Propagate forces from the end to the base
% We can model constraint force at the cut loop as external wrench at the
% last body. They are modelled in absolute coordinates (following Phi and
% Phi_q from the HandC algorithm).
f_ext = cell(1,model.NB);

for i = 1:length(model.loop)  
    p = find(model.loop{i}==-1, 1, 'last' ); % cut link 1
    s = find(model.loop{i}==1, 1, 'last' ); % cut link 2
    % apply external wrench to cut link      
    f_ext{p} = f_loop{i};
    if ~isempty(s)
        f_ext{s} = -f_loop{i};
    end
end

f = apply_external_forces( model.parent, Xup, f, f_ext );

% Step 3: Compute total joint force and Transform spatial force to 
% generalized forces at the joints
for i = model.NB:-1:1
  %tau(i,1) = S{i}' * f{i};
  if model.parent(i) ~= 0
    % recurrent joint force Eq.(5.10)  
    % Note: Xup = X(j-1)2(j)   Xup* = X(j-1)2(j)^-T    Xup' = X(j)2(j-1)^T
    f{model.parent(i)} = f{model.parent(i)} + Xup{i}'*f{i};
  end
end

% Compute relative joint forces
f_rel = f;


% Step 4: Transform all forces to the base Plucker coordinates!!
f_abs  = cell(1,model.NB);
for i = 1:model.NB
    if model.parent(i) == 0
      Xa{i} = Xup{i};
    else
      Xa{i} = Xup{i} * Xa{model.parent(i)};
    end
    f_abs{i} =  Xa{i}' * f_rel{i};
end
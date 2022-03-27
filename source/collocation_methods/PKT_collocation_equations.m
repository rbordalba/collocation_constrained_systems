% given a function [xDot, Falgebraic] = f(x,u,z,h), the sampling time h, and 
% the degree d of the collocation polynomial,
% it returns a function Collocation_eq = F.
function F = PKT_collocation_equations(f)

import casadi.*

nx = f.size1_in(0); % number of states
nu = f.size1_in(1); % number of actions
nz = f.size1_in(2); % number of algebraic variables in dyn equation (qddot and lambda)
nl = f.size1_in(3); % number of slack variable for velocity

% create CasADi function for discrete-time dynamics using collocation 
Xk = MX.sym('Xk', nx, 1); % state k
Zk = MX.sym('Zk', nz, 1); % algebraic at k
Xkp = MX.sym('Xkp', nx, 1); % state k+1
Zkp = MX.sym('Zkp', nz, 1); % algebraic at k+1
Zc = MX.sym('Zc', nz, 1); % algebraic at collocation
Uk = MX.sym('Uk', nu, 1); % action at the start of collocation interval
Ukp = MX.sym('Ukp', nu, 1); % action at the end of collocation interval
Gammac = MX.sym('gammac', nl, 1); % extra variable for collocation

% nonlinear equations and cost contribution for collocation interval
eq = {};
h = MX.sym('h',1);

% eval dynamics at knot points
[XkDot, algk] = f(Xk, Uk, Zk, zeros(nl,1));
[XkpDot, algkp] = f(Xkp, Ukp, Zkp, zeros(nl,1));

% collocation and algebraic equations
Uc = (Uk + Ukp)/2;
t = 0.5;
Xc = (2*t^3-3*t^2+1)*Xk + (t^3-2*t^2+t)*h*XkDot + (-2*t^3+3*t^2)*Xkp + (t.^3-t.^2)*h*XkpDot;
[fc, algc] = f(Xc, Uc, Zc, Gammac);

% derivative at collocation 
pcDot = (6*t^2-6*t)*Xk + (3*t^2-4*t+1)*h*XkDot + (-6*t^2+6*t)*Xkp + (3*t^2-2*t)*h*XkpDot;
collocation_eq = h*fc - pcDot;

% append collocation equations
eq = {eq{:}, algk};
eq = {eq{:}, collocation_eq};
eq = {eq{:}, algc};
eq = {eq{:}, algkp};


% implicit discrete-time dynamics
F = Function('F', {Xk, Xkp, Uk, Ukp, Zk, Zkp, Zc, Gammac, h}, ...
                  {vertcat(eq{:})});

end
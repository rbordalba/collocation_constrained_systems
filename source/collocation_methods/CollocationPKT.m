%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
%                           TRAJECTORY SMOOTHING                          % 
%                    (Posa-Kuindersma-Tedrake DirCol)                     %
%                                                                         %
%   Implementation: Ricard Bordalba                                       %
%                                                                         %
%   Institut de Robotica i Informatica Industrial (CSIC-UPC)              %
%   Kinematics and robot design group                                     %
%                                                                         %
%   Description: Trajectory smoothing/optimization with PKT method.       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import casadi.*
% Limits in actuators and joints
p.xmin = p.model.xmin; p.xmin(p.n_q+1:end) = -repmat(p.max_vel,p.n_q,1); 
p.xmax = p.model.xmax; p.xmax(p.n_q+1:end) = repmat(p.max_vel,p.n_q,1);
ts_x_opt.Data(:,p.n_q+1:end) = ts_x_opt.Data(:,p.n_q+1:end); % unwrap angles to have continuous trajectory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              casADi variables                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define differential states
q = SX.sym('q', p.n_q);
qDot = SX.sym('qDot', p.n_q);
x = [q; qDot];
% Lagrange multipliers
lamb = SX.sym('lamb', p.n_l);
% time variable
h = SX.sym('h');
% define controls
u = SX.sym('u', p.n_u);
% slack variable for velocity
gamma = SX.sym('gamma', p.n_l);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              casADi dynamics                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute kinematic and dynamic terms
[M, C, phi, phi_q, xi] = HandC(p.model, q, qDot);

% Velocity constraint equation
phiDot = phi_q*qDot;

% Semi-explicit dynamics (lambda is an algebraic variable)
qDDot = M\(p.model.B*u-p.model.friction.*qDot-C-phi_q'*lamb);

% differential equation
xDot = [qDot+phi_q'*gamma; qDDot];

% kinematic constraints
Fkin = [phi;
        phiDot;
        phi_q*qDDot-xi];
    
% fwd dynamics Function
xDotFunction = Function('xDotFunction', {x, u, lamb, gamma}, {xDot});

% kinematic function at knot points
kinematic_functions = Function('kinematic_functions', {[q; qDot],u, lamb}, {Fkin});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           casADi collocation                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collocation variables
Xk = SX.sym('Xk', p.n_x, 1); % state k
Lambdak = SX.sym('Lambdak', p.n_l, 1); % lagrange multipliers
Xkp = SX.sym('Xkp', p.n_x, 1); % state k+1
Lambdakp = SX.sym('Lambdakp', p.n_l, 1); % lagrange multipliers
Lambdac = SX.sym('Lambdac', p.n_l, 1); % lagrange multipliers at collocation
Gammac = SX.sym('gammac', p.n_l, 1); % slack variable for collocation
Uk = SX.sym('Uk', p.n_u, 1); % action at the start of collocation interval
Ukp = SX.sym('Ukp', p.n_u, 1); % action at the end of collocation interval

% eval dynamics at knot points
XkDot = xDotFunction(Xk, Uk, Lambdak, zeros(p.n_l,1));
XkpDot = xDotFunction(Xkp, Ukp, Lambdakp, zeros(p.n_l,1));

% collocation state and action
Xc = (Xk + Xkp)/2 + h*(XkDot - XkpDot)/8;
Uc = (Uk + Ukp)/2;

% dynamics at col
XcDot = xDotFunction(Xc, Uc, Lambdac, Gammac);
% polynomial derivative at col
pcDot = -3*(Xk - Xkp)/2 - (XkDot + XkpDot)*h/4;
% collocation
collocationEq = h*XcDot - pcDot;
collocationFunction = Function('collocationFunction', {Xk, Xkp, Uk, Ukp, Lambdak, Lambdakp, Lambdac, Gammac, h}, {collocationEq});               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            casADi NLP formulation                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];
P={};
P0 = [];

% Fixed initial state
Xk = xs;

% Action at 0
Uk = MX.sym('U_0', p.n_u);
w = {w{:}, Uk};
lbw = [lbw; p.umin];
ubw = [ubw; p.umax];
w0 = [w0;  getdatasamples(ts_u_opt,1)']; % use initial guess

Lambdak = MX.sym('Lambda_0', p.n_l, 1);
w = {w{:}, Lambdak};
lbw = [lbw; -inf(p.n_l,1)];
ubw = [ubw; inf(p.n_l,1)];
% explicitly compute lambda from initial guess
xc = getdatasamples(ts_x_opt,1)';
uc = getdatasamples(ts_u_opt,1)';
[~, ~, lambda0] = ForwardDynamics([],xc,uc,p.model,3);
w0 = [w0; lambda0];

% manifold constraints
kin_eq = kinematic_functions(Xk,Uk,Lambdak);
g = {g{:}, kin_eq(end-p.n_l+1:end)}; % only acceleration cnstr
lbg = [lbg; zeros(p.n_l,1)];
ubg = [ubg; zeros(p.n_l,1)]; 

% Formulate the NLP
for k=1:N
    
    % New NLP variable for the time
    hK = MX.sym(['h_' num2str(k-1)]);
    w = {w{:}, hK};
    if opt.minimize_time || opt.free_time
        lbw = [lbw; 0.1*h0(k)]; % time as variable then
        ubw = [ubw; 1*h0(k)];    
    else
        lbw = [lbw; 1*h0(k)];
        ubw = [ubw; 1*h0(k)];
    end
    w0 = [w0; h0(k)];
    if k>1
        g = {g{:}, hK-hprev};
        lbg = [lbg; 0];
        ubg = [ubg; 0]; 
    end
    
    % New NLP variable for the algebraic variables in the collocation
    Lambdac = MX.sym(['Lambdac_' num2str(k)], p.n_l, 1);
    w = {w{:}, Lambdac};
    lbw = [lbw; -inf*ones(p.n_l,1)];
    ubw = [ubw; inf*ones(p.n_l,1)];
    % explicitly compute lambda from initial guess
    xc = getdatasamples(ts_x_opt,k)';
    uc = getdatasamples(ts_u_opt,k)';
    [~, ~, lambda0] = ForwardDynamics([],xc,uc,p.model,3);
    w0 = [w0; lambda0];

    % New NLP variable for PKT collocaiton
    Gammac = MX.sym(['Gammac_' num2str(k)], p.n_l); 
    w = {w{:}, Gammac};
    lbw = [lbw; -opt.max_gamma*ones(p.n_l,1)]; % wierd cost trajectories when increasing gamma
    ubw = [ubw; opt.max_gamma*ones(p.n_l,1)]; % 2
    w0 = [w0; zeros(p.n_l,1)];

    % New NLP variable for state at end of interval
    Xkp = MX.sym(['X_' num2str(k+1)], p.n_x); 
    w = {w{:}, Xkp};
    lbw = [lbw; p.xmin];
    ubw = [ubw;  p.xmax];
    w0 = [w0; getdatasamples(ts_x_opt,k+1)'];
    
    % New NLP variable for the control at k+1
    Ukp = MX.sym(['U_' num2str(k+1)], p.n_u);
    w = {w{:}, Ukp};
    lbw = [lbw;  p.umin];
    ubw = [ubw;  p.umax];
    w0 = [w0;  getdatasamples(ts_u_opt,k+1)']; % action from CUIK trajectory for initial guess
    
    % New variable for algebraic 
    Lambdakp = MX.sym(['Lambda_' num2str(k+1)], p.n_l, 1);
    w = {w{:}, Lambdakp};
    lbw = [lbw; -inf(p.n_l,1)];
    ubw = [ubw; inf(p.n_l,1)];
    % explicitly compute lambda from initial guess
    xc = getdatasamples(ts_x_opt,k+1)';
    uc = getdatasamples(ts_u_opt,k+1)';
    [~, ~, lambda0] = ForwardDynamics([],xc,uc,p.model,3);
    w0 = [w0; lambda0];
    
    % Collocation constraint
    coll_eq = collocationFunction(Xk, Xkp, Uk, Ukp, Lambdak, Lambdakp, Lambdac, Gammac, hK);

    % Collocation
    g = {g{:}, coll_eq};
    lbg = [lbg; zeros(size(coll_eq))];
    ubg = [ubg; zeros(size(coll_eq))]; 

    % manifold constraints
    kin_eq = kinematic_functions(Xkp,Ukp,Lambdakp);
    g = {g{:}, kin_eq};
    lbg = [lbg; zeros(size(kin_eq))];
    ubg = [ubg; zeros(size(kin_eq))]; 
   
    % Minimize u
    % exact integral = h/3 · [ (uk)^2 + (uk+1)^2 + uk · uk+1 ]
    J = J + opt.minimize_u*hK/3*(Uk'*Uk+Ukp'*Ukp+Uk'*Ukp);
    
    % Minimize udot
    J = J + opt.minimize_udot*1/2*(Ukp-Uk)'*(Ukp-Uk)*hK;        
    
    % Minimize time
    J = J + opt.minimize_time*hK;
    
    Xk = Xkp;
    Uk = Ukp;
    Lambdak = Lambdakp;
    hprev = hK;
end

% Terminal constraint
g = [g, {(Xk-xg)}]; % set fixed goal
lbg = [lbg; -1e-8*ones(p.n_x,1)];
ubg = [ubg; 1e-8*ones(p.n_x,1)];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            casADi NLP solver                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
ipopt_opts = struct('linear_solver', 'ma27','tol',opt.tolerance,'constr_viol_tol',opt.constraint_tolerance);
opts = struct('ipopt', ipopt_opts, 'expand', true, 'jit', false);
solver = nlpsol('solver', 'ipopt', prob, opts);
% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, ...
           'lbg', lbg, 'ubg', ubg);
       
w_opt = [xs; full(sol.x)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Solution                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the solution
x_opt = [];
xdot_opt = [];
u_opt = [];
h_opt = [];
t_opt = [];
t_col = [];
t_knot = [];

% unwrap data from solution vector
x_index = 1:p.n_x;
u_index = p.n_x+1:p.n_x+p.n_u;
qddot_index = p.n_x+p.n_u+p.n_l+1:p.n_x+p.n_u+p.n_l+p.n_q;
h_index = p.n_x+p.n_u+p.n_l+1;
jump = p.n_x+p.n_u+p.n_l*3+1;

              
for k=0:N
    % Evaluate only knot points
    xopt = w_opt(x_index+jump*k);
    x_opt = [x_opt xopt];   
        
    % in zoh, N action point, in foh, N+1 action points
    if k==0
        u_opt = [u_opt w_opt(u_index+jump*k)];
    end
    
    if k<N
        
        % N+1 knot points, N intervals
        h_opt = [h_opt w_opt(h_index+jump*k)];
        
        % in zoh, N action point, in foh, N+1 action points
        u_opt = [u_opt w_opt(u_index+jump*(k+1))];
          
        % time at knot points 
        t_knot = [t_knot sum(h_opt(1:k))];        
        t_opt = [t_opt t_knot(end)]; 
        t_col = [t_col t_knot(end) t_knot(end)+h_opt(k+1)/2];

    else
        % time at last knot point
        t_knot = [t_knot sum(h_opt(1:k))];
        t_opt = [t_opt t_knot(end)];
        t_col = [t_col t_knot(end)];
    end
end

% Optimized trajectory to timeseries
ts_x_opt = timeseries(x_opt',t_opt);
ts_x_opt.Data(:,p.n_q+1:end) = ts_x_opt.Data(:,p.n_q+1:end);
ts_u_opt.Time = t_knot; % time at knot points
ts_u_opt.Data = u_opt';

% Store error for each d
stats = solver.stats;
success = stats.success;
exec_time = stats.t_wall_total;
cost = full(sol.f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
%                           TRAJECTORY SMOOTHING                          % 
%               (Naive, Baumgarte or Orthogonal Projection)               %
%                                                                         %
%   Author: Ricard Bordalba                                               %
%                                                                         %
%   Institut de Robotica i Informatica Industrial (CSIC-UPC)              %
%   Kinematics and robot design group                                     %
%                                                                         %
%   Description: Trajectory smoothing/optimization with Basic, Baumgarte, %
%   or orthogonal projection.                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import casadi.*
% Read limits in joints from model file
p.xmin = p.model.xmin; p.xmin(p.n_q+1:end) = -repmat(p.max_vel/p.scaleVel,p.n_q,1); 
p.xmax = p.model.xmax; p.xmax(p.n_q+1:end) = repmat(p.max_vel/p.scaleVel,p.n_q,1);
ts_x_opt.Data(:,p.n_q+1:end) = ts_x_opt.Data(:,p.n_q+1:end)/p.scaleVel; % unwrap angles to have continuous trajectory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              casADi variables                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define differential states
q = SX.sym('q', p.n_q);
qDot = SX.sym('qDot', p.n_q);
x = [q; qDot];
% define algebraic states
lamb = SX.sym('lamb', p.n_l);
% time variable
h = SX.sym('h');
% define controls
u = SX.sym('u', p.n_u);
% qDDot is treated as an algebraic in DAE
qDDot = SX.sym('qDDot', p.n_q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              casADi dynamics                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute kinematic and dynamic terms
[M, C, phi, phi_q, xi, phi_q_dot] = HandC(p.model, q, qDot*p.scaleVel);

% Velocity constraint equation
phiDot = phi_q*qDot;
% Baumgarte stabilization?
if opt.Baumgarte
    alpha = 2; % 1 < alpha < 20 (tunned following [Blajer-CMAME2011-Methods for constraint violation])
    beta = 0.5; 
    xi = xi-alpha*phiDot*p.scaleVel-beta*phi;
end

% Implicit dynamics
z = [lamb; qDDot];
n_z = p.n_l+p.n_q; % number of algebraics in DAE (lambda + qDDot)
dynamics = M*qDDot+C+phi_q'*lamb-p.model.B*u+p.model.friction.*qDot*p.scaleVel; % explicit dynamics (better for DAE)   

% differential equation
xDot = [qDot*p.scaleVel; qDDot/p.scaleVel];

% algebraic equation
Falgebraic = phi_q*qDDot-xi;
Falgebraic = [Falgebraic;
              dynamics];

% jacobian of consistency condition 
JacobianEq = [phi_q zeros(p.n_l, p.n_q); phi_q_dot/p.scaleVel phi_q];

% manifold function - consistency conditions
xManifoldF = Function('xManifoldF', {[q; qDot]}, {[phi; phiDot]});
% jacobian function
xJacobianF = Function('xJacobianF', {[q; qDot]}, {JacobianEq});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            casADi collocation                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continuous time dynamics DAE
f = Function('f', {x, u, z}, {xDot, Falgebraic});
% Collocation function for DAE
F = DAE_foh_collocation_equations(f, opt.d, opt.scheme);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            casADi NLP formulation                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_root = [0 collocation_points(opt.d, opt.scheme)]; % Get collocation points

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

% Formulate the NLP
for k=1:N
    
    % New NLP variable for the time
    if opt.minimize_time || opt.free_time
        hK = MX.sym(['h_' num2str(k-1)]);
        w = {w{:}, hK};
        % free time step if we minimize time
        lbw = [lbw; 0.1*h0(k)];
        ubw = [ubw; 1.5*h0(k)];   
        w0 = [w0; h0(k)];
        if k>1
            g = {g{:}, hK-hprev};
            lbg = [lbg; 0];
            ubg = [ubg; 0]; 
        end    
    else
        hK = h0(1);
    end
    
    % New NLP variable for the state variables in the collocation 
    CXk = MX.sym(['CX_' num2str(k)], p.n_x, opt.d);
    w = [w, {reshape(CXk, opt.d*p.n_x, 1)}];
    lbw = [lbw; repmat(p.xmin,opt.d,1)];
    ubw = [ubw; repmat(p.xmax,opt.d,1)];
    xk0 = getdatasamples(ts_x_opt,k)'; % initial knot state k
    xkp0 = getdatasamples(ts_x_opt,k)'; % initial knot state k+1
    xcol0 = getdatasamples(resample(ts_x_opt,h0(1)*(k-1) + tau_root(2:end)*h0(1)),1:opt.d)';
    for j=1:opt.d
       w0 = [w0; xcol0(:,j)]; 
    end
    
    % New NLP variable for the algebraic variables in the collocation
    CZk = MX.sym(['CZ_' num2str(k)], n_z, opt.d);
    w = [w, {reshape(CZk, opt.d*n_z, 1)}];
    lbw = [lbw; -inf(n_z*opt.d,1)];
    ubw = [ubw; inf(n_z*opt.d,1)];
    % explicitly compute lambda from initial guess
    uk0 = getdatasamples(ts_u_opt,k)'; 
    ukp0 = getdatasamples(ts_u_opt,k+1)'; 
    ucol0 = uk0 + (ukp0-uk0).*tau_root(2:end);
    for j=1:opt.d
        [qDDot0, ~, lambda0] = ForwardDynamics([],xcol0(:,j),ucol0(:,j),p.model,3);
        w0 = [w0; lambda0; qDDot0];
    end
    
    % New NLP variable in case of orthogonal projection to manifold
    if opt.projection
        MUk = MX.sym(['MU_' num2str(k)], p.n_e); 
        w = [w, {MUk}];
        lbw = [lbw; -ones(p.n_e,1)];
        ubw = [ubw; ones(p.n_e,1)];
        w0 = [w0; zeros(p.n_e,1)];
    end

    % New NLP variable for state at end of interval
    Xkp = MX.sym(['X_' num2str(k+1)], p.n_x); 
    w = [w, {Xkp}];
    lbw = [lbw; p.xmin];
    ubw = [ubw;  p.xmax];
    w0 = [w0; getdatasamples(ts_x_opt,k+1)'];
    
    % New NLP variable for the control at k+1
    Ukp = MX.sym(['U_' num2str(k+1)], p.n_u);
    w = {w{:}, Ukp};
    lbw = [lbw;  p.umin];
    ubw = [ubw;  p.umax];
    w0 = [w0;  getdatasamples(ts_u_opt,k+1)']; % action from CUIK trajectory for initial guess
    
    % Collocation constraint
    [coll_eq, Xnext] = F(Xk, reshape(CXk, opt.d*p.n_x, 1), Uk, Ukp, reshape(CZk, n_z*opt.d,1), hK);
    
    % Collocation
    g = {g{:}, coll_eq};
    lbg = [lbg; zeros(n_z*opt.d+p.n_x*opt.d,1)];
    ubg = [ubg; zeros(n_z*opt.d+p.n_x*opt.d,1)]; 

    % Continuity conditions 
    % 1: orthogonal projection, else: no projection
    DF = xJacobianF(Xkp);
    if opt.projection  
        jump = Xkp-Xnext;

        % penalize wierd projections
        J = J + opt.projection_penalty*hK*(jump'*jump);
        
        % using orthogonal correction with extra variable mu        
        continuity = [xManifoldF(Xkp);
                      Xkp-Xnext+DF'*MUk];
    else
        % without using correction
        continuity = Xkp-Xnext;
    end

    g = [g, {continuity}];
    lbg = [lbg; zeros(size(continuity))];
    ubg = [ubg; zeros(size(continuity))];

    % Minimize u
    % exact integral if zoh = h/3 · [ (uk)^2 + (uk+1)^2 + uk · uk+1 ]
    J = J + opt.minimize_u*hK/3*(Uk'*Uk+Ukp'*Ukp+Uk'*Ukp);
    
    % Minimize udot
    J = J + opt.minimize_udot*1/2*(Ukp-Uk)'*(Ukp-Uk)*hK;        
    
    % Minimize time
    J = J + opt.minimize_time*hK;
       
    Xk = Xkp;
    Uk = Ukp;
    if opt.minimize_time || opt.free_time
        hprev = hK;
    end
end

% Terminal constraint
if strcmp(optMode,"Basic") || strcmp(optMode,"Baumgarte")
    if opt.nullspace_goal
        Z = null(full(xJacobianF(xg)));
        g = [g, {Z'*(Xk-xg)}]; % set fixed goal
        lbg = [lbg; -1e-8*ones((p.n_q-p.n_l)*2,1)];
        ubg = [ubg; 1e-8*ones((p.n_q-p.n_l)*2,1)];
    else
        g = [g, {(Xk-xg)}]; % set fixed goal with lower threshold
        lbg = [lbg; -1e-2*ones(p.n_x,1)];
        ubg = [ubg; 1e-2*ones(p.n_x,1)];        
    end     
else
     g = [g, {(Xk-xg)}]; % set fixed goal
     lbg = [lbg; -1e-8*ones(p.n_x,1)];
     ubg = [ubg; 1e-8*ones(p.n_x,1)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Solution                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_opt = [];
u_opt = [];
h_opt = [];
t_opt = [];
t_col = [];
x_col = [];
t_knot = [];
x_knot = [];
% unwrap data from solution vector
x_index = 1:p.n_x;
u_index = p.n_x+1:p.n_x+p.n_u;
if opt.minimize_time || opt.free_time
    h_index = p.n_x+p.n_u+1;
end
xcol_index = p.n_x+p.n_u+1:p.n_x+p.n_u+p.n_x*opt.d;
jump = p.n_x+p.n_u+p.n_x*opt.d+n_z*opt.d;
if opt.projection % extra variable for projection
    jump = jump+p.n_e;
end
if opt.minimize_time || opt.free_time
    jump = jump+1;
    xcol_index = xcol_index + 1;
end

% Recover collocation polynomials
tau_root = [0 collocation_points(opt.d, opt.scheme)]; % Get collocation points
% trick: with Radau, col. and knot point coincide and timeseries reorders data
if strcmp(opt.scheme,'radau')
   tau_root(end) = tau_root(end)-1e-12; 
end
              
for k=0:N
    % Evaluate only knot points
    xopt = w_opt(x_index+jump*k);
    x_opt = [x_opt xopt];   
    
    t_knot = [t_knot sum(h_opt(1:k))];
    x_knot = [x_knot xopt];
    
    % in zoh, N action point, in foh, N+1 action points
    if k==0
        u_opt = [u_opt w_opt(u_index+jump*k)];
    end
    
    if k<N
        
        % N+1 knot points, N intervals
        if opt.minimize_time || opt.free_time
            h_opt = [h_opt w_opt(h_index+jump*k)];
        else
           h_opt = [h_opt h0(1)];
        end
         % in zoh, N action point, in foh, N+1 action points
        u_opt = [u_opt w_opt(u_index+jump*(k+1))];
          
        xCol = reshape(w_opt(xcol_index+jump*k),p.n_x,opt.d); % recover collocation points
        x_opt = [x_opt xCol]; 
        % time at knot points and collocation
        t_opt = [t_opt sum(h_opt(1:k))+tau_root*h_opt(k+1)]; 
        
        t_col = [t_col sum(h_opt(1:k))+tau_root(2:end)*h_opt(k+1)];
        x_col = [x_col xCol];

    else
        % time at last knot point
        t_opt = [t_opt sum(h_opt)];
    end
    
end

% Optimized trajectory to timeseries
ts_x_opt = timeseries(x_opt',t_opt);
ts_x_opt.Data(:,p.n_q+1:end) = ts_x_opt.Data(:,p.n_q+1:end)*p.scaleVel;
knotTime = [0 cumsum(h_opt)];
ts_u_opt.Time = knotTime; % time at knot points
ts_u_opt.Data = u_opt';

% Store error for each d
stats = solver.stats;
success = stats.success;
exec_time = stats.t_wall_total;
cost = full(sol.f);
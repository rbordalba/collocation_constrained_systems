%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
%                           TRAJECTORY SMOOTHING                          % 
%                            (Local coordinates)                          %
%                                                                         %
%   Author: Ricard Bordalba                                               %
%                                                                         %
%   Institut de Robotica i Informatica Industrial (CSIC-UPC)              %
%   Kinematics and robot design group                                     %
%                                                                         %
%   Description: Trajectory smoothing with casADi with local coordinates  %
%                methods. we use foh for action, Legendre polynomials     %
%                to approximate dynamics, and tangent space               %
%                local coordinates to guarantee driftless solutions.      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import casadi.*

y0 = zeros(p.n_y,1);
p.ymin = -p.max_y*ones(p.n_y,1);
p.ymax = p.max_y*ones(p.n_y,1);
% Read limits in joints from model file
p.xmin = p.model.xmin; p.xmin(p.n_q+1:end) = -repmat(p.max_vel/p.scaleVel,p.n_q,1); 
p.xmax = p.model.xmax; p.xmax(p.n_q+1:end) = repmat(p.max_vel/p.scaleVel,p.n_q,1);
ts_x_opt.Data(:,p.n_q+1:end) = ts_x_opt.Data(:,p.n_q+1:end)/p.scaleVel; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              casADi variables                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define differential states
q = SX.sym('q', p.n_q);
qDot = SX.sym('qDot', p.n_q);
x = [q; qDot];
% define tangent space parameters
x_chart = SX.sym('x_chart', p.n_x); % chart/tg space center 
U_chart = SX.sym('U_chart', p.n_x, p.n_y); % chart/tg space basis
y = SX.sym('y', p.n_y);
% time variable
h = SX.sym('h');
% define controls
u = SX.sym('u', p.n_u);
% algebraic variables of implicit dynamics
lamb = SX.sym('lamb', p.n_l);
qDDot = SX.sym('qDDot', p.n_q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              casADi dynamics                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute kinematic and dynamic terms
[M, C, phi, phi_q, xi, phi_q_dot] = HandC(p.model, q, qDot*p.scaleVel);
% constraint equation
phiDot = phi_q*qDot;
% differential equation (implicit form)
dynamics = M*qDDot+C+phi_q'*lamb-p.model.B*u+p.model.friction.*qDot*p.scaleVel;   
% algebraic variables
p.n_z = p.n_x+p.n_l+p.n_q;
z = [x; lamb; qDDot];
% fwd dynamics
xDot = [qDot*p.scaleVel; qDDot/p.scaleVel];
% Dynamics in tg space   
yDot = U_chart'*xDot;
% Implicit inverse map using tg space coordinates
psi = [y - U_chart'*(x-x_chart);
       StateSpace_Manifold(x, p.model)];
% algebraic equation
Falgebraic = [dynamics;
              phi_q*qDDot-xi;
              psi];
% Precompute continuos basis Uc using Gram-Smith orthonormalization
[~, DF] = StateSpace_Manifold(x_chart,p.model);
% Gram Smith orthogonalization 
[~,Unew] = orthonormalize(DF,U_chart);          
          
% CasADi function: Inverse map function
PsiF = Function('PsiF', {y, x, x_chart, U_chart}, {psi});
% CasADi function: orthonormalization function
orthonormalizeF = Function('orthonormalizeF', {x_chart,U_chart}, {Unew});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            casADi collocation                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Continuous time dynamics in local with the corresponding algebraic eq.
f = Function('f', {y, u, z, x_chart, U_chart}, {yDot, Falgebraic});

% Collocation function
F = tangent_foh_collocation_equations(f, opt.d, opt.scheme);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            casADi NLP formulation                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get collocation points
tau_root = [0 collocation_points(opt.d, opt.scheme)]; 

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
p0 = [];

% New NLP variable for initial state
Xk = xs;
Yk = zeros(p.n_y,1);

% Action at 0
Uk = MX.sym('U_0', p.n_u);
w = {w{:}, Uk};
lbw = [lbw; p.umin];
ubw = [ubw; p.umax]; 
w0 = [w0;  getdatasamples(ts_u_opt,1)']; % action from CUIK trajectory for initial guess
    
% Compute Tg space basis
xchart = getdatasamples(ts_x_opt,1)'; % Xk
[~, DF] = StateSpace_Manifold(xchart,p.model);  
Uchart = null(DF);
UchartK = Uchart;
XchartK = xchart;
% Formulate the NLP
for k=1:N

    % New NLP variable for the time step
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
    
    % New NLP variable for the algebraic variables in the collocation 
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
    
    % New NLP variable for the collocation states 
    CYk = MX.sym(['CY_' num2str(k)], p.n_y, opt.d);
    w = [w, {reshape(CYk, opt.d*p.n_y, 1)}];
    lbw = [lbw; repmat(p.ymin,opt.d,1)];
    ubw = [ubw; repmat(p.ymax,opt.d,1)];
    for j=1:opt.d
       w0 = [w0; Uchart'*(xcol0(:,j)-xk0)]; % mid col. local coordinates between k and k+1
    end
    
    CLk = MX.sym(['CL_' num2str(k)], p.n_l+p.n_q, opt.d);
    w = [w, {reshape(CLk, opt.d*(p.n_l+p.n_q), 1)}];
    lbw = [lbw; -inf((p.n_l+p.n_q)*opt.d,1)];
    ubw = [ubw; inf((p.n_l+p.n_q)*opt.d,1)];      
    % explicitly compute lambda from initial guess
    xc = getdatasamples(ts_x_opt,k)';
    ucol0 = getdatasamples(resample(ts_u_opt,h0(1)*(k-1) + tau_root(2:end)*h0(1)),1:opt.d)';
    uc = getdatasamples(ts_u_opt,k)'; 
    [qDDot0, ~, lambda0] = ForwardDynamics([],xc,uc,p.model,3);
    lambda0 = [lambda0; qDDot0];
    w0 = [w0; repmat(lambda0,opt.d,1)];    
    % append algebraic
    CZk = [CXk; CLk];
    
    % New NLP variable for local coordinates at k+1 in chart k
    Ykp = MX.sym(['Y_' num2str(k+1)], p.n_y);
    w = {w{:}, Ykp};
    lbw = [lbw;  p.ymin];
    ubw = [ubw;  p.ymax];
    w0 = [w0; Uchart'*(getdatasamples(ts_x_opt,k+1)'-xc)]; 
    
    % New NLP variable for state at end of interval
    Xkp = MX.sym(['X_' num2str(k+1)], p.n_x); % x(k+1)
    w = [w, {Xkp}];
    lbw = [lbw; p.xmin];
    ubw = [ubw;  p.xmax];
    xkp = getdatasamples(ts_x_opt,k+1)'; % estimated x(k+1)
    w0 = [w0; xkp];
    
    if opt.variable_chart
        % start new polynomial at chart centre Xchartk = Xk (so Yk = 0)
        Yk = zeros(p.n_y,1); 
    else
        % when chart is fixed, start of polynomial can vary
        Yk = UchartK'*(Xk-XchartK); 
    end
    % New NLP variable for the control at k+1
    Ukp = MX.sym(['U_' num2str(k+1)], p.n_u);
    w = {w{:}, Ukp};
    lbw = [lbw;  p.umin];
    ubw = [ubw;  p.umax];
    w0 = [w0;  getdatasamples(ts_u_opt,k+1)']; 
    
    % Collocation constraint
    [coll_eq, Yk_end] = F(Yk, reshape(CYk, opt.d*p.n_y, 1), Uk, Ukp, ...
    reshape(CZk, p.n_z*opt.d,1), XchartK, UchartK, hK);
              
    % Collocation conditions
    g = {g{:}, coll_eq};
    lbg = [lbg; zeros((p.n_y+p.n_z)*opt.d,1)];
    ubg = [ubg; zeros((p.n_y+p.n_z)*opt.d,1)];
    
    % Continuity conditions: x(k+1) = psi_k(y(end)) using chart k
    g = [g, {Ykp-Yk_end}];
    lbg = [lbg; zeros(p.n_y,1)];
    ubg = [ubg; zeros(p.n_y,1)];
    
    if k==N
        % manifold constraint at the last point is implicitly defined when
        % xN = xg
        Gn = PsiF(Ykp, Xkp, XchartK, UchartK);
        g = [g, {Gn(1:p.n_y)}];
        lbg = [lbg; zeros(p.n_y,1)];
        ubg = [ubg; zeros(p.n_y,1)];
    else
        g = [g, {PsiF(Ykp, Xkp, XchartK, UchartK)}];
        lbg = [lbg; zeros(p.n_x,1)];
        ubg = [ubg; zeros(p.n_x,1)];
    end
    
    % avoid wierd projections
    er = XchartK+UchartK*Ykp-Xkp;
    J = J + opt.local_penalty*hK*(er'*er);
    if opt.extra_avoid_proj
        jump = er'*er;
        g = {g{:}, jump};
        lbg = [lbg; 0];
        ubg = [ubg; 0.5];
    end
    
    % Compute Tg space basis
    xchart = getdatasamples(ts_x_opt,k+1)'; % Xk
    [~, DF] = StateSpace_Manifold(xchart,p.model);
    Uchart = null(DF); % basis of nullspace
    if opt.variable_chart
        UchartKp = orthonormalizeF(Xkp,Uchart);
        XchartKp = Xkp;
    else
        UchartKp = Uchart; 
        XchartKp = xchart;
    end
    
    % Cost function
    % Minimize u: exact integral = h/3 · [ (uk)^2 + (uk+1)^2 + uk · uk+1 ]
    J = J + opt.minimize_u*hK/3*(Uk'*Uk+Ukp'*Ukp+Uk'*Ukp);
    
    % Minimize udot
    J = J + opt.minimize_udot*1/2*(Ukp-Uk)'*(Ukp-Uk)*hK;        
    
    % Minimize time
    J = J + opt.minimize_time*hK;
        
    Xk = Xkp;
    Uk = Ukp;
    UchartK = UchartKp;
    XchartK = XchartKp;
    if opt.minimize_time || opt.free_time
        hprev = hK;
    end
end

% Terminal constraint 
g = [g, {(Xk-xg)}];
lbg = [lbg; -1e-8*ones(p.n_x,1)];
ubg = [ubg; 1e-8*ones(p.n_x,1)];
    
% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}), 'p', vertcat(P{:}));
ipopt_opts = struct('linear_solver', 'ma27','tol',opt.tolerance,'constr_viol_tol',opt.constraint_tolerance);
opts = struct('ipopt', ipopt_opts, 'expand', true, 'jit', false);
solver = nlpsol('solver', 'ipopt', prob, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     Results                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
costf = Function('costf', {vertcat(w{:}),vertcat(P{:})},{J});
constraintf = Function('constraintf', {vertcat(w{:}),vertcat(P{:})},{vertcat(g{:})});
 
% compute initial cost
cost = full(costf(w0,p0));

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'p', p0, ...
            'lbg', lbg, 'ubg', ubg);
w_opt = [xs;full(sol.x)];

%% decode ipopt solution
x_opt = [];
y_opt = [];
u_opt = [];
h_opt = [];
t_opt = [];
t_col = []; 
Uchart_opt = zeros(p.n_x,p.n_y,N+1);
xchart_opt = zeros(p.n_x,N+1);

% index info
x_index = 1:p.n_x;
u_index = p.n_x+1:p.n_x+p.n_u;
if opt.minimize_time || opt.free_time
    h_index = p.n_x+p.n_u+1;
end
xcol_index = p.n_x+p.n_u+1:p.n_x+p.n_u+opt.d*p.n_x;
ycol_index = p.n_x+p.n_u+opt.d*p.n_x+1:p.n_x+p.n_u+opt.d*p.n_x+opt.d*p.n_y;
jump = p.n_y*(opt.d+1)+p.n_u+opt.d*p.n_z+p.n_x;
if opt.minimize_time || opt.free_time
    jump = jump+1;
    xcol_index = xcol_index + 1;
    ycol_index = ycol_index + 1;
end

% Recover collocation polynomials
tau_root = [0 collocation_points(opt.d, opt.scheme)]; % Get collocation points
% trick: with Radau, col. and knot point coincide and timeseries reorders data
if strcmp(opt.scheme,'radau')
   tau_root(end) = tau_root(end)-1e-14; 
end
costi = 0;
for k=0:N

    % state knot point
    xopt = w_opt(x_index+jump*k);
    x_opt = [x_opt xopt];
    
    % Compute Tg space basis from initial guess
    xchart = getdatasamples(ts_x_opt,k+1)'; % Xk
    [~, DF] = StateSpace_Manifold(xchart,p.model);
    Uchart = null(DF);
    if k~=0 && opt.variable_chart
        xchart_opt(:,k+1) = xopt; % store charts
        % Chart orthonormalization with respect to initial guess
        Uchart_opt(:,:,k+1) = full(orthonormalizeF(xchart_opt(:,k+1),Uchart));
    else
        xchart_opt(:,k+1) = xchart; % store charts
        Uchart_opt(:,:,k+1) = Uchart;
    end
    
    % paremeter knot point
    if opt.variable_chart
        yopt = zeros(p.n_y,1); 
    else
        yopt = Uchart'*(xopt-xchart); % fixed value
    end
    y_opt = [y_opt yopt];

    % in zoh, N action point, in foh, N+1 action points
    if k==0
        u_opt = [u_opt w_opt(u_index+jump*k)];
    end       
    

    % collocation points
    if k<N
        % N+1 knot points, N intervals
        if opt.minimize_time || opt.free_time
            h_opt = [h_opt w_opt(h_index+jump*k)];
        else
            h_opt = [h_opt h0(1)];
        end
         % in zoh, N action point, in foh, N+1 action points
        u_opt = [u_opt w_opt(u_index+jump*(k+1))];
        costi = costi+h_opt(end)/3*(u_opt(:,end)'*u_opt(:,end)+u_opt(:,end)'*u_opt(:,end-1)+u_opt(:,end-1)'*u_opt(:,end-1));
        
        % Collocation points
        xCol = reshape(w_opt(xcol_index+jump*k),p.n_x,opt.d); % recover collocation points
        x_opt = [x_opt xCol];      
        yCol = reshape(w_opt(ycol_index+jump*k),p.n_y,opt.d); % recover collocation points
        y_opt = [y_opt yCol];

        % time at knot points and collocation
        t_opt = [t_opt sum(h_opt(1:k))]; 
        t_opt = [t_opt t_opt(end)+tau_root(2:end)*h_opt(k+1)]; 
        
        t_col = [t_col sum(h_opt(1:k))+tau_root(2:end)*h_opt(k+1)];

    else
        % time at last knot point
        t_opt = [t_opt sum(h_opt)];
    end
    

end

%cost = [cost costi]; % actually measure cost related to u'u
cost = [cost full(sol.f)];

% Use optimized trajectory to compute new charts in next iteration
ts_x_opt = timeseries(x_opt', t_opt);
ts_x_opt.Data(:,p.n_q+1:end) = ts_x_opt.Data(:,p.n_q+1:end)*p.scaleVel;
ts_y_opt = timeseries(y_opt', t_opt);
h0 = h_opt; 
knotTime = [0 cumsum(h_opt)];
ts_u_opt.Time = knotTime; % time at knot points,
ts_u_opt.Data = u_opt';

ts_Uchart_opt = timeseries(Uchart_opt,knotTime); 
ts_xchart_opt = timeseries(xchart_opt',knotTime); % not scaled!

% Store error for each d
stats = solver.stats;
success = stats.success;
exec_time = stats.t_wall_total;
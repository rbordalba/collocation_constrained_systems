%ForwardDynamics   Computes the acceleration using system dynamics.
%   qdd = ForwardDynamics(t,x,u,model,integration,f_ext) returns the acceleration. 
%
%   OUTPUTS:
%       - qdd System acceleration.
%   INPUTS:
%       - t: time of simulation. In general t=[], since dynamics will not
%       depend on time t. Its needed when using ode solvers such as ode45.
%       - x: state (q,qdot).
%       - u: input forces vector. IMPORTANT:
%       In constrained systems, size(u) is not equal to size(q) and use the
%       projection matrix B to obtain tau=Bu.
%       In open chain robots, tau = u.
%       - model: model function from Featherstone Toolbox.
%       - integration: determines the method to compute the forward
%       dynamics of the robot in constrained systems! Not needed for open
%       chain robots.
%           1. Baumgarte stabilization method for the constraints.
%           2. Constraint Violation Supression (CVSU) method.
%           3. Simple ODE
%           4. Coordiante partitioning using actuated joints.
%           5. Projection method from Aghili
%           6. Coordiante partitioning using Tg. space local coordinates
%       - f_ext: spatial or planar force vector. See Featherstone toolbox
%       notes for more info about it.
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu

function [qdd, Fc, lambda] = ForwardDynamics(t,x,u,model,integration,f_ext)
n_q = length(x)/2;
n_u = model.n_u;
n_e = model.n_e;
B = model.B;

q = x(1:n_q); %position
v = x(n_q+1:2*n_q); %velocity

if  isa(u, 'function_handle')
   u = u(t,x);
end
% Convert u to proper dimensions
u = B*u;
% Viscous and Coulomb friction
if isfield(model,'friction')
    u = u-model.friction.*v;
end
if isfield(model,'Coulomb')
    vlim = 0.01;
    u = u-model.Coulomb.*tanh(v/vlim); %sigmoid Coulomb
end

% Constrained system 
if isfield(model,'loop')  
    % Openloop dynamic terms and Constraints equations
    if nargin == 6 
        if isa(f_ext, 'function_handle')
            f_ext=f_ext(t,x);
        end
        [M, f, Phi, Phi_q, xi] = HandC(model,q,v,f_ext);
    else
        [M, f, Phi, Phi_q, xi] = HandC(model,q,v);
    end
    
    f = u-f;
    n_lambda = length(xi);
   
    % Forward dynamics for constrained systems
    switch integration
        
        case 1 % Baumgarte
            tau = 0.15;
            alpha = 2/tau;
            Phi_dot = Phi_q*v;
            xi = xi -alpha*Phi_dot-alpha^2/4*Phi;
            temp = ([Phi_q' M; zeros(n_lambda) Phi_q]\[f; xi]); % pinv?
            lambda = temp(1:n_lambda);
            Fc = Phi_q'*lambda;
            qdd = temp(n_lambda+1:end);
        case 2 %CVSU method for constraint violation suppression
            Ts = 0.01;
            Phi_dot = Phi_q*v;
            temp = ([Phi_q' M; zeros(n_lambda) Phi_q]\[f; xi-Phi_dot/Ts]);
            lambda = temp(1:n_lambda);
            Fc = Phi_q'*lambda;
            qdd=temp(n_lambda+1:end);
        case 3 % Common extended mass matrix approach
            temp = ([Phi_q' M; zeros(n_lambda) Phi_q]\[f; xi]);
            lambda = temp(1:n_lambda);
            Fc = Phi_q'*lambda;
            qdd = temp(n_lambda+1:end);
        case 4 % ODE with coordinate partitioning
            % Choosing actuated joints as independent
            Qi = B;
            dep = find(~any(B, 2));   
            Qd = eye(n_q); Qd = Qd(:,dep);
            Phi_qi = Phi_q*Qi; 
            Phi_qd = Phi_q*Qd;
            A = Qi*eye(n_u);
            A(dep,:)=-Phi_qd\Phi_qi;
            d = -Phi_q\xi;
            ddqi = (A'*M*A)\(A'*(f+M*d));
            qdd = A*ddqi-d;
        case 5 % ODE with projection method from Aghili standard
            P = eye(size(Phi_q,2))-pinv(Phi_q)*Phi_q;
            Cq_dot = -pinv(Phi_q)*-xi;
            M_tilde = P*M-(P*M)'; %skew symmetric
            Mc = M+M_tilde;
            qdd = (Mc)\(P*f+M*Cq_dot);
            Fc = (eye(size(P))-P)*(f-M*qdd);   
            lambda = []; 
        case 6 %coordiante partitioning using tangent space coordinates(C-space Haug)
            n_y=model.NB-n_e/2;
            G_q = [Phi_q; null(Phi_q)'];
            G_y = [zeros(n_e/2,n_y);-eye(n_y)];
            D = -G_q\G_y; %D=null(Phi_q) directly if we use x as center of chart
            Gqdotqd = [-xi;zeros(n_y,1)];
            d = (G_q\Gqdotqd);
            ydd = (D'*M*D)\(D'*(f+M*d));
            qdd = D*ydd-d;
    end
    
% Open-chain system
else
    qdd = FDab(model,q,v,u);
end
end
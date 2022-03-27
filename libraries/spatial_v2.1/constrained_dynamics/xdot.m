 %XDOT   Computing state transition dynamics 
%   dxdt = XDOT(t,x,u,model,integration) returns the state dynamics.
%
%   OUTPUTS:
%       - dxdt state derivative x_dot = (qdot,qdotdot)
%   INPUTS:
%       - t: time of simulation. In general t=[], since dynamics will not
%       depend on time t.
%       - x: state.
%       - u: input forces vector. IMPORTANT:
%       In constrained systems, size(u) is not equal to size(q), we use the
%       projection matrix B to obtain tau=Bu.
%       - model: model function from Featherstone Toolbox.
%       - integration: determines the method to compute the forward
%       dynamics of the robot (only when model is a constrained system).
%           1. Baumgarte stabilization method for the constraints.
%           2. Simple ODE
%           3. Constraint Violation Supression (CVSU) method.
%           4. Coordiante partitioning using actuated joints.
%           5. Projection method from Aghili
%           6. Coordiante partitioning using Tg. space local coordinates
%       - f_ext: spatial or planar force vector. See Featherstone toolbox
%       notes for more info about it.
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu

function dxdt = xdot(t,x,u,model,integration,f_ext)
n_q = length(x)/2;
q = x(1:n_q); %position
qd = x(n_q+1:2*n_q); %velocity

if nargin == 6
    qdd = ForwardDynamics(t,x,u,model,integration,f_ext);
else
    qdd = ForwardDynamics(t,x,u,model,integration);
end

dxdt = [qd;qdd];
end


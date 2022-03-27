%YDOT   Computing state transition dynamics in parameter space.
%   dydt = YDOT(t,y,xc,u,model,integration) returns the state dynamics in
%   parameter space (tangent space parametrization with chart center at
%   xc). Note that it can fail if the parameter y is far away from the
%   validity of the parametrization. This means that we can compute x from
%   the given y at the corresponding chart.
%
%   OUTPUTS:
%       - dydt state derivative in parameter space.
%   INPUTS:
%       - t: time of simulation. In general t=[], since dynamics will not
%       depend on time t.
%       - y: state in parameter space.
%       - xc: center of the chart of the tangent space parametrization.
%       - u: input forces vector.
%       - model: model function from Featherstone Toolbox.
%       - integration: determines the method to compute the forward
%       dynamics of the robot.
%           1. Simple ODE
%           2. Projection method from Aghili
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu

function dydt = ydot(t,y,xc,u,model,integration)

    if  isa(u, 'function_handle')
       u = u(t,y);
    end

    % Project to manifold
    [x, ~] = InverseMap(y,xc.x0,xc.x0,xc.U0,model);
    
    dxdt = xdot(t,x,u,model,integration);
    
    dydt = xc.U0'*dxdt;
    
end


%INVERSEMAP   Inverse map from parameter to ambient speed in a
%             constrained system using tangent space parametrization.
%   x = InverseVelocityMap(yDot, xc, Uc, model) returns the ambient speed xdot
%   corresponding to the parameter speed ydot from the chart whose center 
%   is xc and its basis is Uc.
%
%   OUTPUTS:
%       - xdot: state.
%   INPUTS:
%       - yDot: state speed in tangent space parametrization.
%       - x: current state.
%       - Uc: tangent space basis at xc. 
%       - model: model function from Featherstone Toolbox.
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu

function xDot = InverseAccelerationMap(yDDot, x, Uc, model)   

    n_y = length(yDot);
    [~, ~, phi, phi_q, xi, phi_q_dot] = HandC(p.model, q, qDot*scaleVel);

    [F, DF] = StateSpace_Manifold(x, model);
    n_e = length(F);

    xDot = [DF; Uc']\[zeros(n_e,n_y); eye(n_y)]*yDot;
end

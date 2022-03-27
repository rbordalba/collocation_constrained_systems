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

function xDot = InverseVelocityMapCspace(yDot, x, Uc, model)   
    
    n_y = length(yDot);
    n_x = length(x);
    q = x(1:n_x/2);
    qDot = x(n_x/2+1:end);
    
    [~,~,Phi,Phi_q,xi,Phi_q_dot] = HandC( model, q, qDot);
    Gq = [Phi_q; Uc'];
    Gy = [zeros(length(Phi),n_y/2); -eye(n_y/2)];
    D = -Gq\Gy;

    xDot(1:n_x/2,1) = D*yDot(1:n_y/2);
    xDot(n_x/2+1:n_x,1) = D*yDot(n_y/2+1:end)-Gq\[Phi_q_dot;zeros(size(Uc'))]*xDot(1:n_x/2,1);
    
end

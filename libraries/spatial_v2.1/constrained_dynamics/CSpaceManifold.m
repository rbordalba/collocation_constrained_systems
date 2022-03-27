%StateSpace_Manifold   Computes the state space manifold X error and its
%Jacobian.
%   [F, DF] = StateSpace_Manifold(x,model) returns F(x) and DF(x). 
%
%   OUTPUTS:
%       - F State Space error F=[Phi;Phi_dot].
%       - DF State Space Jacobian
%   INPUTS:
%       - x: state (q,qdot).
%       - model: model function from Featherstone Toolbox.
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu
function [F, DF] = CSpaceManifold(q,model)
    n_q = length(q);
    
    if isfield(model,'loop')
        [~,~,Phi,Phi_q] = HandC(model,q,zeros(n_q,1));
    else
       warning('Model should include a loop field to describe the closed-chain.');
       return;
    end

    F = Phi;
    DF = Phi_q;
      
end
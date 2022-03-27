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
function [F, DF] = StateSpace_Manifold(x,model)
    n_q = length(x)/2;
    q = x(1:n_q);
    v = x(n_q+1:end);
    
    if isfield(model,'loop')
        [~,~,Phi,Phi_q,xi,Phi_q_dot] = HandC(model,q,v);
    else
       warning('Model should include a loop field to describe the closed-chain.');
       return;
    end
    
    Phi_dot = Phi_q*v;     
    F = [Phi;Phi_dot];

    % Compute the contrained Jacobian of F = [Phi; Phi_dot]
    DF = [Phi_q     zeros(size(Phi_q));
          Phi_q_dot             Phi_q];   

      
      
end
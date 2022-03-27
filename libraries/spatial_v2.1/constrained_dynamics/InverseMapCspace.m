%INVERSEMAP   Inverse map from parameter to ambient coordiantes in a
%             constrained system using tangent space parametrization.
%   x = InverseMap(y,x0,xc,U0,model) returns the ambient coordinates x
%   corresponding to the parameter y from the chart whose center is xc and
%   its basis is U0.
%
%   OUTPUTS:
%       - x: state.
%   INPUTS:
%       - y: state parameter in tangent space parametrization.
%       - xc: chart center.
%       - U0: tangent space basis at xc. 
%       - model: model function from Featherstone Toolbox.
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu

function [x, exitflag] = InverseMapCspace(y, x0, xc, Uc, model)   
nx = length(x0);
nq = nx/2;
ny = length(y);
qc = xc(1:nq);
q = x0(1:nq);
% find q implicitly
[Phi, Phi_q] = CSpaceManifold([q;zeros(nq,1)], model);
H =  [Phi;
      Uc'*(q-qc)-y(1:ny/2)];     

iter = 0;  exitflag = 1;  
while norm(H)>1e-08
    % Newton step
    B = [Phi_q;Uc'];
    s = -B\H;
    q = q + s; 

    % Check new value
    [Phi, Phi_q] = CSpaceManifold([q;zeros(nq,1)], model);
    H = [Phi;
         Uc'*(q-qc)-y(1:ny/2)];  
    
    iter = iter+1;
    if iter > 30
        q = x0(1:nq);%zeros(size(x0));
        exitflag = 0; % failed to converge
        %warning('Newton diverged');
        break;
    end
end
%%
if exitflag <= 0
    Phi = @(q) CSpaceManifold([q;zeros(nq,1)], model);
    H =  @(q)[Phi(q);
              Uc'*(q-qc)-y(1:ny/2)];     
    options = optimoptions('fsolve','Display','none');
    q0 = x0(1:nq);
    [q, exitflag] = fsolve(H,q0,options);
    if exitflag <= 0
        q = q0;
    end
end

% find q dot explicitlly
[Phi, Phi_q] = CSpaceManifold([q;zeros(nq,1)], model);
Gq = [Phi_q; Uc'];
Gy = [zeros(length(Phi),ny/2); -eye(ny/2)];
D = -Gq\Gy;
qDot = D*y(ny/2+1:end);

x = [q;qDot];

end
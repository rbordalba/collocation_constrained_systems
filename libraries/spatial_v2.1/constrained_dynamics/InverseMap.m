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

function [x, exitflag] = InverseMap(y, x0, xc, Uc, model)   
x = x0; % Initial guess

% Simple Newton method
if false

    [F, DF] = StateSpace_Manifold(x, model);

    H = [F;
         Uc'*(x-xc)-y];     

    iter = 0;  exitflag = 1;  
    while norm(H)>1e-15
        % Newton step
        B = [DF;Uc'];
        s = -B\H;
        x = x + s; 

        % Check new value
        [F, DF] = StateSpace_Manifold(x, model);
        H = [F;
             Uc'*(x-xc)-y];  

        iter = iter+1;
        if iter > 10
            x = x0;
            exitflag = 0; % failed to converge
            %warning('Newton diverged');
            break;
        end
    end
else
    % Matlab solver
    F = @(x) StateSpace_Manifold(x, model);
    H =  @(x)[F(x);
              Uc'*(x-xc)-y];     
    options = optimoptions('fsolve','Display','none','OptimalityTolerance',1e-15);
    [x,~, exitflag] = fsolve(H,x0,options);
    if exitflag <= 0
        x = x0;
    end
end



end
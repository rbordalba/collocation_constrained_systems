% given a function [xDot, Falgebraic] = f(x,u,z,h), the sampling time h, and 
% the degree d of the collocation polynomial, and the collocation scheme,
% it returns a function [Collocation_eq, xnext] = F.
function F = DAE_foh_collocation_equations(f, d, scheme)

import casadi.*

nx = f.size1_in(0);
nz = f.size1_out(1);
nu = f.size1_in(1);

% set collocation points
tau_rt = [0 casadi.collocation_points(d, scheme)];   %vector length d+1
tau_root = num2cell(tau_rt);

% coefficients of the collocation equation
C = zeros(d+1,d+1);

% coefficients of the continuity equation
D = zeros(d+1, 1);

% construct polynomial basis
for j=1:d+1
  % construct Lagrange polynomials to get the polynomial basis at the collocation point
  coeff = 1;
  for r=1:d+1
    if r ~= j     %~= : not equal
      coeff = conv(coeff, [1, -tau_root{r}]);
      coeff = coeff / (tau_root{j}-tau_root{r});
    end
  end
  
  % evaluate the polynomial at the final time to get the coefficients of the continuity equation
  D(j) = polyval(coeff, 1.0);  %evaluates polynom coeff at point 1.0

  % evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  pder = polyder(coeff);     %returns derivative of polynomial coeff
  for r=1:d+1
    C(j,r) = polyval(pder, tau_root{r});
  end
end

% create CasADi function for discrete-time dynamics using collocation 
% (it should speed up jitting. Adapted from 
% https://gist.github.com/jaeandersson/9bf6773414f30539daa32d9c7deeb441)
for j=1:d
    Xc{j} = MX.sym(['Xc_' num2str(j)], nx, 1); % states at collocation
    Zc{j} = MX.sym(['Zc_' num2str(j)], nz, 1); % algebraic at collocation
end
Uk = MX.sym('Uk', nu, 1); % action at the start of collocation interval
Ukp = MX.sym('Ukp', nu, 1); % action at the end of collocation interval

X0_c = MX.sym('X0_c', nx, 1); % state at the beginning of the collocation

% state at the end of the collocation interval
Xf = D(1)*X0_c;
for j=1:d    
    % add contribution to the end state
    Xf = Xf + D(j+1)*Xc{j};
end

% nonlinear equations and cost contribution for collocation interval
eq = {};
h = MX.sym('h',1);

% collocation and algebraic equations
for j = 1:d
    % expression for the state derivative at the collocation point
    xp = C(1,j+1)*X0_c;
    for r=1:d
        xp = xp + C(r+1,j+1)*Xc{r};  % equation 5 ex sheet
    end
    
    % evaluate the function at the collocation points
    Uc = Uk + tau_rt(j+1)*(Ukp-Uk);
    [fj, algj] = f(Xc{j}, Uc, Zc{j});

    % append collocation equations
    eq = {eq{:}, h*fj - xp};
    eq = {eq{:}, algj};
    
end

% implicit discrete-time dynamics
F = Function('F', {X0_c, vertcat(Xc{:}), Uk, Ukp, vertcat(Zc{:}), h}, ...
                  {vertcat(eq{:}), Xf});

end
% This function creates the collocation equations in the tangent space 
% INPUTS:
%        - A function [yDot, Falgebraic] = f(y, u, z, x_chart, U_chart).
%        where y is the local coord, u is the input, z are the algebraic
%        variables, x_chart is the chart centre and U_chart is the chart
%        basis.
%        It outputs the dynamics in local coord, i.e. yDot, and the
%        algebraic equations required to solve yDot. It includes the
%        implicit map psi(x,y,x_chart,U_chart)=0 and the impl. dynamics.
%        - the degree d of the collocation polynomial.
%        - the collocation scheme, either radau or legendre
% OUTPUT:
%        - A function to create the collocation equations for each interval

function F = tangent_foh_collocation_equations(f, d, scheme)

import casadi.collocation_points
import casadi.MX
import casadi.Function

ny = f.size1_in(0); % dimension of tg coordinates y
nz = f.size1_out(1); % dimension of algebraic variables
nu = f.size1_in(1); % dimension of action coordinates
nx = f.size1_in(3); % dimension of state x

% set collocation points
tau_rt = [0 collocation_points(d, scheme)];   %vector length d+1
tau_root = num2cell(tau_rt);

% coefficients of the collocation equation
C = zeros(d+1,d+1);

% coefficients of the continuity equation
D = zeros(d+1, 1);

% coefficients of the quadrature function
B = zeros(d+1, 1);

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

  % evaluate the integral of the polynomial to get the coefficients of the quadrature function
  pint = polyint(coeff);   % integral of polynom (Integrationskonstante 0)
  B(j) = polyval(pint, 1.0);
  
end

% create CasADi function for discrete-time dynamics using collocation 
% (it should speed up jitting. Adapted from 
% https://gist.github.com/jaeandersson/9bf6773414f30539daa32d9c7deeb441)
for j=1:d
    Yc{j} = MX.sym(['Yc_' num2str(j)], ny, 1);
    Zc{j} = MX.sym(['Zc_' num2str(j)], nz, 1);
end
Uk = MX.sym('Uk', nu, 1); % action at the start of collocation interval
Ukp = MX.sym('Ukp', nu, 1); % action at the end of collocation interval
Y0_c = MX.sym('Y0_c', ny, 1); % state at the beginning of the collocation

% state at the end of the collocation interval
Yf = D(1)*Y0_c;
for j=1:d    
    % add contribution to the end state
    Yf = Yf + D(j+1)*Yc{j};
   
% % % %     % Add contribution to quadrature function
% % % %     J = J + B(j+1)*qj*dt;
       
end

% nonlinear equations and cost contribution for collocation interval
eq = {};
Qf = 0;
Rf = [];
h = MX.sym('h',1);
x_chart = MX.sym('x_chart', nx); % chart/tg space center 
U_chart = MX.sym('U_chart', nx, ny); % chart/tg space basis

% collocation and algebraic equations
for j = 1:d
    % expression for the state derivative at the collocation point
    yp = C(1,j+1)*Y0_c;
    for r=1:d
        yp = yp + C(r+1,j+1)*Yc{r};  % equation 5 ex sheet
    end
    
    % evaluate the function at the collocation points
    Uc = Uk + tau_rt(j+1)*(Ukp-Uk);
    [fj, algj] = f(Yc{j}, Uc, Zc{j}, x_chart, U_chart);

    % append collocation equations
    eq = {eq{:}, h*fj - yp};
    eq = {eq{:}, algj};
end

% implicit discrete-time dynamics
F = Function('F', {Y0_c, vertcat(Yc{:}), Uk, Ukp, vertcat(Zc{:}), x_chart, U_chart, h}, ...
                  {vertcat(eq{:}), Yf});
end
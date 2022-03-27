%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Function to reconstruct state trajectory with barycentric interpolation %
% INPUTS:                                                                 %
%    t:    Time to evaluate                                               %
%    ts_x: Timeseries that should include knot points and collocation     %
%          points x(1) ... x(k) x(k,1) ... x(k,d) x(k+1) ... x(N)         %
%    d:    Legendre order. Order of interpolated polynomial for the       %
%          intervals is d+1.                                              %
% scheme:  Collocation scheme, either radau or legendre.                  %
%   out:  The default output is the state x. Set out=2 if we want xDot.   %
%                                                                         %
% OUTPUT:                                                                 %
%   x:    State at time t                                                 %
%   xDot: State derivative at time t                                      %
%                                                                         %
%   Author: Ricard Bordalba                                               %
%                                                                         %
%   Institut de Robotica i Informatica Industrial (CSIC-UPC)              %
%   Kinematics and robot design group                                     %
%   2020                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, xDot] = stateTrajectory(t, ts_x, d, scheme, outputSelect)
% Recover collocation polynomials (range [0-1])
tau_root = [0 casadi.collocation_points(d, scheme)]; % Get collocation points

% Find the k polynomial corresponding to time t.
%tPos = find(abs(ts_x.Time-t)<1e-12); 
tPos = find(ts_x.Time == t);
%tPos = find(ts_x.Time == t);
if isempty(tPos) 
    tmp = sort([t; ts_x.Time]);
    tPos = find(tmp==t);
    k = floor((tPos-2)/(d+1)); % interval of selectred time (starting at 0)
    knotPoint = [];
else % case its a knot point
    tPos = min(tPos);
    k = floor((tPos-1)/(d+1));
    knotPoint = tPos-k*(d+1);
    % In case it's the end of trajectory, we use spline N-1
    if length(ts_x.Time) == tPos && strcmp(scheme,'legendre')
        k = k - 1;
        knotPoint = [];
    end    
end

% take spline states and times
vk = ts_x.Data((1:d+1)+k*(d+1),:)'; % states defining k poly
tk = ts_x.Time((1:d+1)+k*(d+1)); % time tk
tkp = ts_x.Time(d+1+k*(d+1)+1); % time tk+1
tEval = (t-tk(1));
At = (tkp-tk(1)); 
% scale time to [0-At]
tau_root = tau_root*At;

% Barycentric interpolation
    
for j=1:d+1  
  % Construct Lagrange polynomials to get the polynomial basis at the collocation point
  coeff(j) = 1/prod(tau_root(j)-tau_root(1:end~=j));  % 1/prod_k~j(x_j-x_k)
end
    
%Compute D
D = zeros(d+1);
for i=1:d+1
    for j=1:d+1
        if i~=j
           D(i,j) = coeff(j)/coeff(i)/(tau_root(i)-tau_root(j));
        end
    end
end
% diagonal case
D(1:d+2:end) = -sum(D,2); 
% differentiate vk
vkDot = (D*vk')';
% differentiate vkDot
vkDDot = (D*vkDot')';

% evaluate with barycentric formula
% case its a knot point
if isempty(knotPoint) 
    den = sum(coeff./(tEval-tau_root));
    num = sum((coeff.*vk)./(tEval-tau_root),2);
    x(:,1) = num./den;
    num = sum(coeff.*vkDot./(tEval-tau_root),2);
    xDot(:,1) = num./den;
    num = sum(coeff.*vkDDot./(tEval-tau_root),2);
    xDDot(:,1) = num./den;
else
   x = vk(:,knotPoint); 
   xDot = vkDot(:,knotPoint); 
   xDDot = vkDDot(:,knotPoint); 
end
    

if nargin == 5 % output select 
    if outputSelect == 2
        x = xDot;
    elseif outputSelect == 3
        x = xDDot;
    else
        x = x;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Function to reconstruct local trajectory with barycentric interpolation %
% INPUTS:                                                                 %
%   t:     Time to evaluate                                               %
%   ts_y:  Timeseries that should include knot points and collocation     %
%          points in local y(k) y(k,1) ... y(k,d) y(k+1) y(k+1,1)         %
%   ts_x:  Timeseries that should include knot points and collocation     %
%          points in ambient x(k) x(k,1) ... x(k,d) x(k+1) x(k+1,1)       %
%   ts_Uc: Timeseries that should include knot points and collocation     %
%          points x(k) x(k,1) ... x(k,d) x(k+1) x(k+1,1)                  %
%   d:     Legendre order. Order of interpolated polynomial for the       %
%          intervals is d+1.                                              %
%  scheme: Collocation scheme, either radau or legendre.                  %
%                                                                         %
% OUTPUT:                                                                 %
%   y:    Local coordinates at time t                                     %
%                                                                         %
%   Author: Ricard Bordalba                                               %
%                                                                         %
%   Institut de Robotica i Informatica Industrial (CSIC-UPC)              %
%   Kinematics and robot design group                                     %
%   2020                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y, yDot, yDDot, xchart, Uchart, x0] = localTrajectory(t, ts_y, ts_x, ts_x_chart, ts_Uc, d, scheme, outputSel)
if t>ts_y.Time(end) || t<ts_y.Time(1)
    error('Time out of bounds');
end

% use knot points as initial guess
x0 = getdatasamples(resample(ts_x,t),1)';

% Recover collocation polynomials (range 0-1)
tau_root = [0 casadi.collocation_points(d, scheme )]; % Get collocation points

tPos = find(abs(ts_y.Time-t)<1e-12); %find(ts_x.Time == t);
if isempty(tPos)
    tmp = sort([t; ts_y.Time]);
    tPos = find(tmp==t);
    k = floor((tPos-2)/(d+1));
    knotPoint = [];
else % knot point
    tPos = min(tPos);
    k = floor((tPos-1)/(d+1)); % interval
    knotPoint = tPos-k*(d+1);
    if length(ts_y.Time) == tPos
        knotPoint = [];
        k = k - 1;
        vk = ts_y.Data(1+k*(d+1),:)';
        tk = ts_y.Time((1:d+1)+k*(d+1));
    end
end



%xchart = ts_x_chart.Data(1+k*(d+1),:)';
xchart = ts_x_chart.Data(k+1,:)';
Uchart = ts_Uc.Data(:,:,k+1);

% take spline points
vk = ts_y.Data((1:d+1)+k*(d+1),:)';
tk = ts_y.Time((1:d+1)+k*(d+1));
tkp = ts_y.Time(d+1+k*(d+1)+1);
tEval = (t-tk(1));
h = (tkp-tk(1));
% scale time
tau_root = tau_root*h;

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
    y(:,1) = num./den;
    num = sum(coeff.*vkDot./(tEval-tau_root),2);
    yDot(:,1) = num./den;
    num = sum(coeff.*vkDDot./(tEval-tau_root),2);
    yDDot(:,1) = num./den;

else
   y = vk(:,knotPoint); 
   yDot = vkDot(:,knotPoint); 
   yDDot = vkDDot(:,knotPoint); 
end

% derivative if outputSel==2
if nargin==8
    if outputSel==2
        y = yDot;
    elseif outputSel==3
        y = yDDot;
    end
end


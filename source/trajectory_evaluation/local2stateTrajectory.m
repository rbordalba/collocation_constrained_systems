%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Function to reconstruct state trajectory from local trajectory          %
% INPUTS:                                                                 %
%   t:     Time to evaluate                                               %
%   ts_y:  Timeseries that should include knot points and collocation     %
%          points in local y(k) y(k,1) ... y(k,d) y(k+1) y(k+1,1)         %
%   ts_x:  Timeseries that should include knot points and collocation     %
%          points in ambient x(k) x(k,1) ... x(k,d) x(k+1) x(k+1,1)       %
%   ts_Uchart: Timeseries that should include knot points and collocation     %
%          points x(k) x(k,1) ... x(k,d) x(k+1) x(k+1,1)                  %
%   d:     Legendre order. Order of interpolated polynomial for the       %
%          intervals is d+1.                                              %
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

function [x, xDot] = local2stateTrajectory(t, ts_y, ts_x, ts_xchart, ts_Uchart, d, scale_vel, scheme, model, outputSel)

[y, yDot, yDDot, xchart, Uchart, x0] = localTrajectory(t, ts_y, ts_x, ts_xchart, ts_Uchart, d, scheme);
    
x = InverseMap(y, x0, xchart, Uchart, model);
    
xDot = InverseVelocityMap(yDot, x, Uchart, model); 

n_q = length(x)/2;
xDot = [ones(n_q,1); scale_vel*ones(n_q,1)].*xDot;
x = [ones(n_q,1); scale_vel*ones(n_q,1)].*x;

if nargin == 10
    if outputSel == 2
        x = xDot;
    end
end
end


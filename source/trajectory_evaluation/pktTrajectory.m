%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Function to reconstruct local trajectory with barycentric interpolation %
% INPUTS:                                                                 %
%   t:        Time to evaluate                                            %
%   ts_x:     Timeseries that should include knot points and collocation  %
%             points in ambient x(k) x(k,1) ... x(k,d) x(k+1) x(k+1,1)    %
%   ts_xdot:  Timeseries that should include knot points and collocation  %
%             points in ambient x(k) x(k,1) ... x(k,d) x(k+1) x(k+1,1)    %
%                                                                         %
% OUTPUT:                                                                 %
%   x:    state                                                           %
%                                                                         %
%   Author: Ricard Bordalba                                               %
%                                                                         %
%   Institut de Robotica i Informatica Industrial (CSIC-UPC)              %
%   Kinematics and robot design group                                     %
%   2020                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p, pdot, pdotdot] = pktTrajectory(t, ts_x, ts_u, model, out)
if t>ts_x.Time(end) || t<ts_x.Time(1)
    error('Time out of bounds');
end

% Find interval
h = diff(ts_x.Time(1:2));
interval = find(ts_x.Time<t, 1, 'last' );
if(isempty(interval))
    interval = 1;
end

xk = getdatasamples(ts_x,interval)';
uk = getdatasamples(ts_u,interval)';
xkdot = xdot(0,xk,uk,model,3);
xkp = getdatasamples(ts_x,interval+1)';
ukp = getdatasamples(ts_u,interval+1)';
xkpdot = xdot(0,xkp,ukp,model,3);

% convert to [0 1] range
t = (t-ts_x.Time(interval))/h;

p = (2*t.^3-3*t.^2+1)*xk + (t.^3-2*t.^2+t)*h*xkdot + ...
    (-2*t.^3+3*t.^2)*xkp + (t.^3-t.^2)*h*xkpdot;

pdot = ((6*t.^2-6*t)*xk + (3*t.^2-4*t+1)*h*xkdot + ...
       (-6*t.^2+6*t)*xkp + (3*t.^2-2*t)*h*xkpdot)/h;
        
pdotdot = ((12*t-6)*xk + (6*t-4)*h*xkdot + ...
          (-12*t+6)*xkp + (6*t-2)*h*xkpdot)/h;
   
% derivative if outputSel==2
if nargin==5
    if out==2
        p = pdot;
    elseif out == 3
        p = pdotdot;
    end
end


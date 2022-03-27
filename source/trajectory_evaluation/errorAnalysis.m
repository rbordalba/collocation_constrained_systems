%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
%                           ERROR ANALYSIS                                % 
%                                                                         %
%   Author: Ricard Bordalba                                               %
%                                                                         %
%   Institut de Robotica i Informatica Industrial (CSIC-UPC)              %
%   Kinematics and robot design group                                     %
%   2020                                                                  %
%                                                                         %
%   Description: script to compute errors in optimization                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

knotTime = t_knot; % knot points
if strcmp(optMode,'Projection')
    knotTimePrime = knotTime(2:end)-1e-12; % knot point prime (for projection only)
else
    knotTimePrime = [];
end
colTime = t_col;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Spline functions                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input
uFunction = @(t) getdatasamples(resample(ts_u_opt,t,'linear'),1)'; %foh

% state trajectory
if strcmp(optMode,'Local')
    ts_x_opt_scaled = ts_x_opt;
    ts_x_opt_scaled.Data(:,p.n_q+1:end) = ts_x_opt_scaled.Data(:,p.n_q+1:end)/p.scaleVel;
    
    yFunction = @(t) localTrajectory(t, ts_y_opt, ts_x_opt_scaled, ts_xchart_opt, ts_Uchart_opt, opt.d, opt.scheme); % y trajectory
    xChartFunction = @(t) getdatasamples(resample(ts_xchart_opt,t,'zoh'),1)'; % chart centre trajectory (constant between intervals)
    UChartFunction = @(t) getdatasamples(resample(ts_Uchart_opt,t,'zoh'),1); % basis of chart (constant between intervals)
    
    xFunction = @(t) local2stateTrajectory(t, ts_y_opt, ts_x_opt_scaled, ts_xchart_opt, ts_Uchart_opt, opt.d, p.scaleVel, opt.scheme, p.model);
    xDotFunction = @(t) local2stateTrajectory(t, ts_y_opt, ts_x_opt_scaled, ts_xchart_opt, ts_Uchart_opt, opt.d, p.scaleVel, opt.scheme, p.model, 2);

elseif strcmp(optMode,'PKT') 
    xFunction = @(t) pktTrajectory(t,ts_x_opt,ts_u_opt,p.model);
    xDotFunction = @(t) pktTrajectory(t,ts_x_opt,ts_u_opt,p.model,2);
    xDotDotFunction = @(t) pktTrajectory(t,ts_x_opt,ts_u_opt,p.model,3);
    qFunction = @(t) eye(p.n_q,p.n_x)*xFunction(t);
    qDotFunction = @(t) eye(p.n_q,p.n_x)*xDotFunction(t);
    qDDotFunction = @(t) eye(p.n_q,p.n_x)*xDDotFunction(t);
    
else
    xFunction = @(t) stateTrajectory(t,ts_x_opt,opt.d,opt.scheme);
    xDotFunction = @(t)stateTrajectory(t,ts_x_opt,opt.d,opt.scheme,2);    
    xDDotFunction = @(t) stateTrajectory(t,ts_x_opt,opt.d,opt.scheme,3);
    qFunction = @(t) eye(p.n_q,p.n_x)*xFunction(t);
    qDotFunction = @(t) eye(p.n_q,p.n_x)*xDotFunction(t);
    qDDotFunction = @(t) eye(p.n_q,p.n_x)*xDDotFunction(t);
end

% dynamic and kinematic errors (ed and ek)
kinematicErrorFunction = @(t) StateSpace_Manifold(xFunction(t), p.model); % ek
dynamicErrorFunction = @(t) xDotFunction(t) - xdot(0, xFunction(t), uFunction(t), p.model, 3); % ed

% Function evaluations at finer time step
hsim = 1e-3; % desired constant time stamp
tsim = linspace(0,ts_x_opt.Time(end),round(ts_x_opt.Time(end)/hsim));
tsim = sort(unique([knotTimePrime knotTime colTime tsim])); % inlcude knot and col points

% evaluate splines
X = [];
XDOT = [];
FWDYN = [];
kinematic_error = [];
dynamic_error = [];
xprev = xs;
for i=1:length(tsim)    
    [xsim, xdotsim] = xFunction(tsim(i));
    usim = uFunction(tsim(i));
    if isnan(xsim)
       pause; 
    end
    X = [X xsim];
    kinematic_error = [kinematic_error StateSpace_Manifold(xsim, p.model)];
    dynamic_error = [dynamic_error xdotsim - xdot(0, xsim, usim, p.model, 3)];
end        

% store finer sampled data
ts_kinematic_error = timeseries(vecnorm(kinematic_error)',tsim);
ts_dynamic_error = timeseries(vecnorm(dynamic_error)',tsim);
ts_x_opt2 = timeseries(X',tsim);
ts_q_opt = timeseries(X(1:p.n_q,:)',tsim);
ts_qDot_opt = timeseries(X(p.n_q+1:end,:)',tsim);

% compute EK and ED from ek and ed
if ~acurate_integrals
    % less accurate integral with trapz
    kinematic_error = trapz(ts_kinematic_error.Time,ts_kinematic_error.Data)/ts_kinematic_error.Time(end);
    dynamic_error = trapz(ts_dynamic_error.Time,ts_dynamic_error.Data)/ts_dynamic_error.Time(end);
else
    % accurate Integral of errors
    kinematic_error = integral(@(t) cellfun(@(t) norm(kinematicErrorFunction(t)),num2cell(t)),0,ts_x_opt.Time(end)-1e-12,'Waypoints',knotTime)/ts_x_opt.Time(end);
    dynamic_error = integral(@(t) cellfun(@(t) norm(dynamicErrorFunction(t)),num2cell(t)),0,ts_x.Time(end)-1e-12,'Waypoints',knotTime)/ts_x_opt.Time(end);
end

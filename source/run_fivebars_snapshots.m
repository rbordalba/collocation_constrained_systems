%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
%                   TRAJECTORY OPTIMIZATION FOR SCHERBOT                  % 
%                                                                         %
%   Author: Ricard Bordalba                                               %
%                                                                         %
%   Institut de Robotica i Informatica Industrial (CSIC-UPC)              %
%   Kinematics and robot design group                                     %
%                                                                         %
%   Required libraries: casADi + IPOPT + MA27 for trajectory opt.         %
%                       spatial_v2 (Featherstone) for dynamic quantities. %
%                                                                         %
%   Description: script to visualize trajectories from  the different     %
%                transcription methods in the fivebars robot              %
%                (aka Scherbot).                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              "Includes"                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paths for required libraries
libraries_path = '../libraries/';
addpath(genpath(libraries_path)); 
addpath(genpath('./')); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Optimization parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optModeList = ["Basic","Baumgarte","PKT","Projection","Local"];
d_list = 3; % list of degree d of the interpolating polynomial
h0 = 0.05; % approx. sampling time for optimization   

opt.scheme = 'legendre'; % collocation scheme: 'radau' or 'legendre'
opt.free_time = false; 
opt.minimize_u = 1; % weight for cost
opt.minimize_udot = 0; % weight for cost
opt.minimize_time = 0; % weight for cost
opt.nullspace_goal = true; % use nullspace goal when there is too much drift 
opt.variable_chart = true; % use orthonormalization to update charts during opt
opt.tolerance = 1e-9; % optimization tolerance
opt.constraint_tolerance = 1e-9; % constraint tolerance in optimization
opt.extra_avoid_proj = false; % ineq constraint to limit projections
opt.projection_penalty = 6e-1; % penalty of the projection in the proj method (larger when d is small and h0 is high)
opt.local_penalty = 9e-1; % penalty of the projection in the local method
opt.max_gamma = 60; % gamma for pkt

% results options
acurate_integrals = false; % uses integral() or an approximation with trapz() to compute errors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Robot parameters                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.model = Scherbot2; % model 
p.max_vel = 40; 
p.max_y = 5;
p.scaleVel = 4; % scaling speed (only optimization side)
p.max_torque = [1.4; 1.4];
p.umin = -p.max_torque; 
p.umax = p.max_torque;
p.n_q = p.model.NB; % number of gen. coordinates (joint angles)
p.n_u = p.model.n_u; % number of actuators
p.n_e = p.model.n_e; % number of pos and vel equations (constraints in X)
p.n_l = p.n_e/2; % number of lagrange multipliers
p.n_x = p.n_q*2; % number of states
p.n_y = 2*p.n_q-p.n_e; % dimension of the tangent space
p.n_p = p.n_u; % number of model parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Load trajectory from CUIK                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert .sol from CUIK to Matlab timeseries
trajectoryPath = './initial_guesses/FiveBars_traj.sol';
[ts_q, ts_dq, ts_u, ts_x] = sol2ts(trajectoryPath, 1000);    
tdesired = ts_x.Time(end);

% Prepare trajectory for initial guess
ts_x.Data(:,1:p.n_q) = unwrap(ts_x.Data(:,1:p.n_q)); % unwrap angles to have continuous trajectory
xs = ts_x.Data(1,:)'; % start state
xg = ts_x.Data(end,:)'; % goal state
t_f = ts_x.Time(end); % final time
N = ceil(t_f/h0); % number of control intervals
h0 = t_f/N; % correct sampling to be exact
topt = 0:h0:t_f; % time at knot points
h0 = diff(topt);
N = length(h0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Trajectory smoothing with casADi                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run for different optimization modes
for optMode = optModeList
        
    opt.d = d_list(1); % Degree of interpolating polynomial

    % resample initial guess at knot points
    ts_x_opt = resample(ts_x, topt);
    ts_u_opt = resample(ts_u, topt);

    % Optimization/ Smoothing mode
    switch optMode
        case "Basic"
            opt.Baumgarte = false;
            opt.projection = false;
            CollocationProjection;
        case "Baumgarte"
            opt.Baumgarte = true;
            opt.projection = false;
            CollocationProjection;
        case "Projection"
            opt.Baumgarte = false;
            opt.projection = true;
            CollocationProjection;
        case "Local"
            CollocationLocal;
        case "PKT"
            CollocationPKT;
        otherwise
            error('Wrong optimization mode')
    end

    % compute errors
    errorAnalysis;

    % store results for plotting
    optModeIndex = find(optModeList==optMode);
    % Store timeseries
    s.ts_x_opt{optModeIndex,opt.d} = ts_x_opt; % knot and col points
    s.ts_u_opt{optModeIndex,opt.d} = ts_u_opt;
    s.ts_x_opt2{optModeIndex,opt.d} = ts_x_opt2; % finer sampling
    s.ts_q_opt{optModeIndex,opt.d} = ts_q_opt; % finer sampling
    s.ts_qDot_opt{optModeIndex,opt.d} = ts_qDot_opt; % finer sampling
    s.ts_kinematic_error{optModeIndex,opt.d} = ts_kinematic_error;
    s.ts_dynamic_error{optModeIndex,opt.d} = ts_dynamic_error;
    
    % Store summary results
    s.stats{optModeIndex,opt.d} = solver.stats;
    s.success(optModeIndex,opt.d) = success;
    s.cost(optModeIndex,opt.d) = cost(end);
    s.execTime(optModeIndex,opt.d) = stats.t_wall_total;
    s.trajTime(optModeIndex,opt.d) = ts_x_opt.Time(end);
    s.dynamic_error(optModeIndex,opt.d) = dynamic_error;
    s.kinematic_error(optModeIndex,opt.d) = kinematic_error;
    s.iterCount(optModeIndex,opt.d) = stats.iter_count;
    s.colTime{optModeIndex,opt.d} = colTime;
    s.knotTime{optModeIndex,opt.d} = knotTime;
    s.knotTimePrime{optModeIndex,opt.d} = knotTimePrime;  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Visualize                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(optModeList)
    videoScherbot(s.ts_x_opt2{i,d_list},0.5,p.model,optModeList(i));
end
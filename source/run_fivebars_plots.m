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
%   Description: script to compare kinematic and dynamic errors for the   %
%                different transcription methods.                         %  
%                                                                         % 
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

% opt options
optModeList = ["Basic","Baumgarte","PKT","Projection","Local"];
opt.d = 3; % degree d of the interpolating polynomial
opt.scheme = 'legendre'; % collocation scheme: 'radau' or 'legendre'
opt.free_time = false; 
opt.minimize_u = 1;
opt.minimize_udot = 0;
opt.minimize_time = 0;
opt.nullspace_goal = true;
opt.variable_chart = true;
opt.h0 = 0.034; % approx. sampling time for optimization   
opt.max_gamma = 60; % gamma for pkt
opt.projection_penalty = 0e-1; % penalization of the projection in the cost (larger when d is small and h0 is high)
opt.local_penalty = 0e-1;
opt.extra_avoid_proj = false;
opt.tolerance = 1e-9;
opt.constraint_tolerance = 1e-9;

% results options
acurate_integrals = false; % accurate integrals with integral() or approximate with trapz()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Robot parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.model = Scherbot2; % model 
p.n_q = p.model.NB; % number of gen. coordinates (joint angles)
p.n_u = p.model.n_u; % number of actuators
p.n_e = p.model.n_e; % number of pos and vel equations (constraints in X)
p.n_l = p.n_e/2; % number of lagrange multipliers
p.n_x = p.n_q*2; % number of states
p.n_y = 2*p.n_q-p.n_e; % dimension of the tangent space
p.n_p = p.n_u; % number of model parameters
p.max_torque = [1.4; 1.4];
p.umin = -p.max_torque; 
p.umax = p.max_torque;
p.max_vel = 40; 
p.max_y = 6;
p.scaleVel = 4; % scaling speed (only optimization side)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Load trajectory from CUIK                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert .sol from CUIK to Matlab timeseries (sampling at 500Hz)
trajectoryPath = './initial_guesses/FiveBars_traj.sol';
[ts_q, ts_dq, ts_u, ts_x] = sol2ts(trajectoryPath, 1000);    
tdesired = ts_x.Time(end);

% Prepare trajectory for initial guess
ts_x.Data(:,1:p.n_q) = unwrap(ts_x.Data(:,1:p.n_q)); % unwrap angles to have continuous trajectory
xs = ts_x.Data(1,:)';
xg = ts_x.Data(end,:)';

t_f = ts_x.Time(end); % final time
N = ceil(t_f/opt.h0); % number of control intervals
h0 = t_f/N; % correct sampling to be exact
topt = 0:h0:t_f;
h0 = diff(topt);
N = length(h0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Trajectory smoothing with casADi                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run for different optimization modes
for optMode = optModeList
        
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
%                              Plotting                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_state_and_action_trajectory = false;
plot_kinematic_error = true;
plot_dynamic_error = true;
plot_convergence = false;
d_plot = opt.d;
run_plots;
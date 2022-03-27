%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
%                   TRAJECTORY OPTIMIZATION FOR SCHERBOT                  % 
%                                                                         %
%   Author: Ricard Bordalba                                               %
%                                                                         %
%   Institut de Robotica i Informatica Industrial (CSIC-UPC)              %
%   Kinematics and robot design group                                     %
%   2020                                                                  %
%                                                                         %
%                                                                         %
%   Required libraries: casADi + IPOPT + MA27 for trajectory opt.         %
%                       spatial_v2 (Featherstone) for dynamic quantities. %
%                                                                         %
%   Description: script to compare the four transcription methods for     %
%                the smoothing problem in the Scherbot.                   %
%                                                                         %
%   Notes: when using projection or local, you may need to tune the       %
%   penalization of the projection step in the CollocationLocal_Scherbot  %
%   and CollocationProjection_Scherbot files.                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Options                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print_latex_table = false;
% paths for required libraries
libraries_path = '../libraries/';
addpath(genpath(libraries_path)); 
% opt options
optModeList = ["Basic","Baumgarte","PKT","Projection","Local"];%,"Baumgarte","PKT","Projection","Local"];
scheme = 'legendre'; % collocation scheme: 'radau' or 'legendre'
free_time = false; 
minimize_u = 1;
minimize_udot = 0;
minimize_time = 0;
nullspace_goal = true;
reduce_vars = false;
variable_chart = true;
d_list = 3; % degree d of the interpolating polynomial
tolerance = 1e-8;
constraint_tolerance = 1e-8;
opt.extraAvoidProj = false; % extra penalize wierd projections in case we have them
opt.h0 = 0.01; % approx. sampling time for optimization   
opt.max_gamma = 60; % gamma for pkt

% results options
acurate_integrals = false; % approximate integrals with integral() or trapz()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Problem parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.model = Scherbot2; % model 

opt.projCost = 0e-1; % penalization of the projection in the cost (larger when d is small and h0 is high)
projection_penalty = 0;
local_penalty = 0e-1;
trajectoryPath = './trajectories/FiveBars_traj2.sol';
p.max_vel = 40; 
p.max_y = 6;%%6.5;
scaleVel = 4; 
p.max_torque = [1.4; 1.4];
p.umin = -p.max_torque; 
p.umax = p.max_torque;

% model and problem dimension
p.n_q = p.model.NB; % number of gen. coordinates (joint angles)
p.n_u = p.model.n_u; % number of actuators
p.n_e = p.model.n_e; % number of pos and vel equations (constraints in X)
p.n_l = p.n_e/2; % number of lagrange multipliers
p.n_x = p.n_q*2; % number of states
p.n_y = 2*p.n_q-p.n_e; % dimension of the tangent space
p.n_p = p.n_u; % number of model parameters


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Load trajectory from CUIK                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert .sol from CUIK to Matlab timeseries (sampling at 500Hz)
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            Baumgarte = false;
            projection = false;
            CollocationProjection;
        case "Baumgarte"
            Baumgarte = true;
            projection = false;
            CollocationProjection;
        case "Projection"
            Baumgarte = false;
            projection = true;
            CollocationProjection;
        case "Local"
            CollocationLocal;
        case "LocalCspace"
            CollocationLocalCspace;
        case "PKT"
            %CollocationPKT;
            CollocationPKTsimplified;
        case "Hermite"
            CollocationHermite;
        otherwise
            error('Wrong optimization mode')
    end

    % store results for plotting
    optModeIndex = find(optModeList==optMode);
    % Store summary results
    s.stats{optModeIndex,opt.d} = solver.stats;
    s.success(optModeIndex,opt.d) = success;
    s.execTime(optModeIndex,opt.d) = stats.t_wall_total;
    s.iterCount(optModeIndex,opt.d) = stats.iter_count;
    s.stats{optModeIndex,opt.d} = solver.stats;

end

%%
plot_state_and_action_trajectory = false;
plot_kinematic_error = false;
plot_dynamic_error = false;
plot_convergence = true;
plot_integration_error = false;
plot_separated_dynamic_error = false;

d_plot = d_list(1);
run_plots;
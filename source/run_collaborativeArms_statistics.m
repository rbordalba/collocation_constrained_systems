%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
%                       TRAJECTORY OPTIMIZATION                           % 
%                                                                         %
%   Author: Ricard Bordalba                                               %
%                                                                         %
%   Institut de Robotica i Informatica Industrial (CSIC-UPC)              %
%   Kinematics and robot design group                                     %
%                                                                         %
%   Required libraries: casADi + IPOPT + MA27 for trajectory opt.         %
%                       spatial_v2 (Featherstone) for dynamic quantities. %
%                                                                         %
%   Description: script to obtain summary table of the different methods  %
%                with different d in the two collaborative arms.          %
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

optModeList= ["Basic","Baumgarte","PKT","Projection","Local"];
d_list = 2:4; % list with degree d of the interpolating polynomial
opt.scheme = 'legendre'; % collocation scheme: 'radau' or 'legendre'
opt.tolerance = 1e-6;
opt.constraint_tolerance = 1e-6;
opt.free_time = false; 
opt.minimize_u = 0;
opt.minimize_udot = 1;
opt.minimize_time = 0;
opt.nullspace_goal = true;
opt.variable_chart = false;
opt.projection_penalty = 1e-1;
opt.local_penalty = 1e-1;
opt.extra_avoid_proj = true;
opt.h0 = 0.04; % approx. sampling time for optimization   
opt.max_gamma = 2; % gamma for pkt

% results options
acurate_integrals = false; % accurate integrals with integral() or approximate with trapz()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.model = fullTwoPandas; % model 
p.max_torque = [87; 87; 87; 87; 12; 12; 12; 87; 87; 87; 87; 12; 12; 12];  
p.max_vel = 8;
p.max_y = 2.8;
p.scaleVel = 1;

% model and problem dimension
p.n_q = p.model.NB; % number of gen. coordinates (joint angles)
p.n_u = p.model.n_u; % number of actuators
p.n_e = p.model.n_e; % number of pos and vel equations (constraints in X)
p.n_l = p.n_e/2; % number of lagrange multipliers
p.n_x = p.n_q*2; % number of states
p.n_y = 2*p.n_q-p.n_e; % dimension of the tangent space
p.n_p = p.n_u; % number of model parameters

% plan param
p.umin = -p.max_torque; 
p.umax = p.max_torque;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Load trajectory from CUIK                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert .sol from CUIK to Matlab timeseries (sampling at 500Hz)
trajectoryPath = './initial_guesses/FullTwoPandas_traj.sol';
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
pkt_done = false;
% run for different optimization modes
for optMode = optModeList
        
    % run for different degree d of polynomial
    for d_index = d_list
        opt.d = d_index; % Degree of interpolating polynomial
    
        % resample at lower rate
        ts_x_opt = resample(ts_x, topt);
        ts_u_opt = resample(ts_u, topt);
        
        % Optimization/ Smoothing mode
        switch optMode
            case "Basic"
                opt.Baumgarte = false;
                opt.projection = false;
                CollocationProjection;
                errorAnalysis;
            case "Baumgarte"
                opt.Baumgarte = true;
                opt.projection = false;
                CollocationProjection;
                errorAnalysis;
            case "Projection"
                opt.Baumgarte = false;
                opt.projection = true;
                CollocationProjection;
                errorAnalysis;
            case "Local"
                % use projection as initial guess for local and resample at opt rate
                ts_x_opt = resample(s.ts_x_opt2{1,d_index},topt);
                ts_u_opt = resample(s.ts_u_opt{1,d_index},topt);
                CollocationLocal;
                errorAnalysis;
            case "PKT"
                % only run PKT once
                if ~pkt_done
                    CollocationPKT;
                    errorAnalysis;
                    pkt_done = true;
                end
            otherwise
                error('Wrong optimization mode')
        end

        if success
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
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Results                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n              Method | Ek·1e-3 |   Ed   |   cost   |   topt     |   iter   \n');

for optMode = optModeList % optmode
    fprintf('--------------------------------------------------------------------------\n');
    optModeIndex = find(optModeList==optMode); 
 
    if strcmp(optMode,'PKT')
        fprintf('%+20s |',num2str(optMode)); 
        fprintf('   %02.2f  |',s.kinematic_error(optModeIndex,d_list(1))*1000); 
        fprintf('  %02.2f  |',s.dynamic_error(optModeIndex,d_list(1))); 
        fprintf('  %5.2f   |',s.cost(optModeIndex,d_list(1))); 
        fprintf('   %4.2f   |',s.execTime(optModeIndex,d_list(1))); 
        fprintf('%5d \n',s.iterCount(optModeIndex,d_list(1))); 
    end    
    for i=1:length(d_list) 
        if ~strcmp(optMode,'PKT')
            fprintf('%+20s |',[num2str(optMode) ' (d = ' num2str(d_list(i)) ')']); 
            fprintf('   %02.2f  |',s.kinematic_error(optModeIndex,d_list(i))*1000); 
            fprintf('  %02.2f  |',s.dynamic_error(optModeIndex,d_list(i))); 
            fprintf('  %5.2f   |',s.cost(optModeIndex,d_list(i))); 
            fprintf('   %4.2f   |',s.execTime(optModeIndex,d_list(i))); 
            fprintf('%5d \n',s.iterCount(optModeIndex,d_list(i))); 
        end
    end
end  
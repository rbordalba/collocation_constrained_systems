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
%   Required libraries: casADi (optimization) + IPOPT + MA27              %
%                       spatial_v2 modified (Featherstone dynamics)       %
%                                                                         %
%   Description: script to compare solutions with different cost function %
%                with the fivebars robot (aka Scherbot) using the         % 
%                projection method.                                       %
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

optMode = ["Projection"]; 
opt.d = 4; % list with degree d of the interpolating polynomial
opt.h0 = 0.015; % approx. sampling time for optimization   
opt.max_gamma = 60; % gamma for pkt
opt.projection_penalty = 6e-1;  % penalization of the projection in the cost (larger when d is small and h0 is high)
opt.tolerance = 1e-9; % optimization tolerance
opt.constraint_tolerance = 1e-9; % constraint tolerance in optimization
opt.nullspace_goal = false; % use nullspace goal when there is too much drift 
opt.free_time = false; 
opt.scheme = 'legendre'; % collocation scheme: 'radau' or 'legendre'

% results options
acurate_integrals = false; % uses integral() or an approximation with trapz() to compute errors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Robot parameters                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model and problem dimension
p.model = Scherbot2; % model 
p.n_q = p.model.NB; % number of gen. coordinates (joint angles)
p.n_u = p.model.n_u; % number of actuators
p.n_e = p.model.n_e; % number of pos and vel equations (constraints in X)
p.n_l = p.n_e/2; % number of lagrange multipliers
p.n_x = p.n_q*2; % number of states
p.n_y = 2*p.n_q-p.n_e; % dimension of the tangent space
p.n_p = p.n_u; % number of model parameters

% plan param
p.max_torque = 1.4;
p.umin = -repmat(p.max_torque,p.n_u,1); 
p.umax = repmat(p.max_torque,p.n_u,1);
p.max_vel = 40; 
p.scaleVel = 4; % downscale velocity to optimizer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Load trajectory from CUIK                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert .sol from CUIK to Matlab timeseries (sampling at 500Hz)
trajectoryPath = './initial_guesses/FiveBars_traj2.sol';
[ts_q, ts_dq, ts_u, ts_x] = sol2ts(trajectoryPath, 500);    
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

% Optimization options for projection method
opt.Baumgarte = false;
opt.projection = true;
        
% run for different degree d of polynomial
for mini=1:3
        switch mini
            case 1
                opt.minimize_u = 1; % minimize integral u'*u
                opt.minimize_udot = 0; % minimize integral udot'*udot
                opt.minimize_time = 0; % minimize time
            case 2
                opt.minimize_u = 0; % minimize integral u'*u
                opt.minimize_udot = 1; % minimize integral udot'*udot
                opt.minimize_time = 0; % minimize time                
            case 3
                opt.minimize_u = 0; % minimize integral u'*u
                opt.minimize_udot = 0; % minimize integral udot'*udot
                opt.minimize_time = 1; % minimize time
        end
        
        % resample at lower rate
        ts_x_opt = resample(ts_x, topt);
        ts_u_opt = resample(ts_u, topt);
        
        CollocationProjection;

        errorAnalysis;

        if success
            % Store data to compare results
            s.ts_x_opt{mini} = ts_x_opt; % knot and col points
            s.ts_x_opt2{mini} = ts_x_opt2; % finer sampling
            s.ts_u_opt{mini} = ts_u_opt;
        end
end

figure; 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.01, 0.3, 0.95]);
a(1) = subplot(4,1,1);
plot(ts_u)
xlabel('t [s]'), ylabel('u(t) [Nm]'), title('Initial guess')
xlim([0 ts_x.Time(end)])
set(gca,'FontSize',14)
a(2) = subplot(4,1,2);
plot(s.ts_u_opt{1})
xlabel('t [s]'), ylabel('u(t) [Nm]'), title('Minimize u')
xlim([0 ts_x.Time(end)])
set(gca,'FontSize',14)
a(3) = subplot(4,1,3);
plot(s.ts_u_opt{2})
xlabel('t [s]'), ylabel('u(t) [Nm]'), title('Minimize udot')
xlim([0 ts_x.Time(end)])
set(gca,'FontSize',14)
a(4) = subplot(4,1,4);
plot(s.ts_u_opt{3})
xlabel('t [s]'), ylabel('u(t) [Nm]'), title('Minimize time')
xlim([0 ts_x.Time(end)])
set(gca,'FontSize',14)

exportgraphics(gcf,'cost.png','Resolution',300)

%videoTrajectory Generate videos of trajectories
%   videoTrajectory Generate videos of trajectories
%
%   OUTPUTS:
%       - 
%   INPUTS:
%       - ts: timeseries to plot.
%       - playbackSpeed: desired playbackSpeed
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu

function videoTrajectory(ts,playbackSpeed,title_name,knotTime)
% Options
fs = 60; % video frame rate
visualize = true; % show or hide recorded figure 
width = 2;

% Set up the movie.
writerObj = VideoWriter('out','MPEG-4'); % Name it.
%writerObj.FileFormat = 'mp4';
writerObj.FrameRate = fs*playbackSpeed; % How many frames per second.
open(writerObj); 

% OPEN FIG
if visualize
    figure;
else
    figure('visible', 'off');
end
tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile;
pos_fig1 = [40 40 800 400];
set(gcf,'Position',pos_fig1)

% figure options
set(gcf,'color','w');
set(gca,'box','on');
set(gca,'FontSize', 34)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read and resample timeseries data at desired time step for visualization
t_s = 1/fs;
n_samples = ceil(ts.Time(end)*fs);

if nargin > 3
    ts_knot = resample(ts,knotTime);
end

h1 = plot(-1,zeros(size(getdatasamples(ts,1))),'LineWidth',width);%,'b','Linewidth',0.5); % trajectory
hold on; 
% plot knot/collocation points
if nargin > 3
    h2 = plot(-1*ones(1,size(ts,2)),zeros(size(getdatasamples(ts_knot,1))),'k.','MarkerSize',width*6,'MarkerFaceColor','k'); % knot times
end
% detect limits
min_y = min(min(ts.Data));
max_y = max(max(ts.Data));
eps = 0.05*(max_y-min_y);
ylimits = [min_y-eps max_y+eps]; % [-1.2 1.2]; %ylimits = [-1.1*min(min(ts.Data)) 1.1*max(max(ts.Data))];
axis([0 ts.Time(end) ylimits]);
xlabel('$t$ [s]','FontSize',12,'Interpreter','latex');
if nargin > 2
    title(title_name,'FontSize',14,'Interpreter','latex');
end
%ylabel('$u_i(t)$ [Nm]','FontSize',14,'Interpreter','latex');
ylabel('$q_i(t)$ [rad]','FontSize',14,'Interpreter','latex');
set(gca,'FontSize',18)
str = [];
for i = 1:7
   str{i} = ['Joint ' num2str(i+7)]; 
end
legend(str,'FontSize',12,'Interpreter','latex');

% update plot
for i=0:n_samples
    tic;

    idx = (ts.Time<=i/fs);
    tmp = getdatasamples(ts,find(idx==1));
    for j=1:length(h1)
        set(h1(j),'XData',ts.Time(idx),'YData',tmp(:,j));
    end
    
    if nargin > 3
        idx = (ts_knot.Time<=i/fs);
        tmp = getdatasamples(ts_knot,find(idx==1));
        for j=1:length(h2)
            set(h2(j),'XData',ts_knot.Time(idx),'YData',tmp(:,j));
        end
    end
    
    tcurr= toc;
    pause(t_s-tcurr); 

    % save frame
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    
end

close(writerObj); % Saves the movie.

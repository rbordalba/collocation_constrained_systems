%video_Scherbot Generate videos of Scherbot following a trajectory
%   video_Scherbot Generate videos of Scherbot following a trajectory
%
%   OUTPUTS:
%       - 
%   INPUTS:
%       - ts_x: timeseries of state trajectory to follow. 
%       - playbackSpeed: desired playbackSpeed
%
%   Ricard Bordalba, IRI-UPC, rbordalba@iri.upc.edu

function videoScherbot(ts_x,playbackSpeed,model,title_name)
% Options
fs = 48; % video frame rate
prevPoints = 24;
trace = 2; % 0 no trace, 1 full load trace with transparency, 2 trace with dots
color = 2; % 1 dark red, 2 gray
actuated = true; % true for red actuated, false for white

figure;
pos_fig = [0.1 0.1 0.4 0.8];
set(gcf, 'Units', 'Normalized', 'OuterPosition', pos_fig);
set(gcf,'color','w');
box off;
hold on; 
axis equal;
ax = gca;
%PLOT WORKSPACE AND INVERSE SINGULARITY LOCI
plotScherbotWorkspace;
Xmin = min(Xcright);
Xmax = max(Xcleft);
Ymin = min(Ycleft);
Ymax = max(Ycleft);
eps = Xmax/10; 
axis([Xmin-eps Xmax+eps Ymin-eps Ymax+eps]);
set(gca,'visible','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read and resample cuik data at desired time step for visualization
n_samples = floor(ts_x.Time(end)*fs);
t_s = 1/fs;
t = linspace(ts_x.Time(1),ts_x.Time(end),n_samples);
ts_x = resample(ts_x,t);
q = ts_x.Data(:,1:5);
qDot = ts_x.Data(:,6:10);
    
% GET ROBOT X AND Y DATA (real)
X = zeros(size(q,1),1);
X = [X X(:,end)-model.L(2)*cos(q(:,1))];
X = [X X(:,end)-model.L(3)*cos(q(:,1)+q(:,2)) ];
X = [X X(:,end)-model.L(4)*cos(q(:,1)+q(:,2)+q(:,3)) ];
X = [X X(:,end)-model.L(5)*cos(q(:,1)+q(:,2)+q(:,3)+q(:,4)) ];

Y = zeros(size(q,1),1);
Y = [Y Y(:,end)-model.L(2)*sin(q(:,1))];
Y = [Y Y(:,end)-model.L(3)*sin(q(:,1)+q(:,2)) ];
Y = [Y Y(:,end)-model.L(4)*sin(q(:,1)+q(:,2)+q(:,3)) ];
Y = [Y Y(:,end)-model.L(5)*sin(q(:,1)+q(:,2)+q(:,3)+q(:,4)) ];

h=[];h1=[];h2=[];h4=[];h3=[];

% plot load trace
if trace == 1
    if color == 1
        h = scatter(X(1,3),Y(1,3),5000,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.6 0 0],'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
    else
        h = scatter(X(1,3),Y(1,3),5000,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
    end  
elseif trace == 2
    h = plot(X(1,3),Y(1,3),'k.','MarkerSize',24);
end
% plot load
if color == 1
    h2 = plot(X(1,3),Y(1,3),'o','MarkerSize',75,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.6 0 0]);
else
    h2 = plot(X(1,3),Y(1,3),'o','MarkerSize',75,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.5 0.5 0.5]);
end
% plot robot
if actuated % red actuated to see dissassembling
    plot(X(1,[1 5]),Y(1,[1 5]),'o','MarkerSize',30,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0],'LineWidth',2); % plot actuated
else
    plot(X(1,[1 5]),Y(1,[1 5]),'o','MarkerSize',30,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','w','LineWidth',2); % plot actuated
end
h3 = plot(X(1,:),Y(1,:),'-','LineWidth',30,'Color',[0 0 0],'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]);
h4 = plot(X(1,:),Y(1,:),'o','MarkerSize',30,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'LineWidth',2); % plot passive

title(title_name,'FontSize',24);
set(findall(gca, 'type', 'text'), 'visible', 'on')

for i=1:size(X,1)
    tic;
    % plot robot trace
    if prevPoints >= i
        tmp = 1:i;
    else
        tmp = i-prevPoints:i;
    end
    set(h,'XData',X(tmp,3),'YData',Y(tmp,3));   
    % plot load
    set(h2,'XData',X(i,3),'YData',Y(i,3));
    % update robot
    set(h3,'XData',X(i,:),'YData',Y(i,:));
    set(h4,'XData',X(i,:),'YData',Y(i,:));
    tcurr= toc;
    pause((t_s-tcurr)/playbackSpeed); 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         % 
%                           SCHERBOT RESULTS                              % 
%                                                                         %
%   Author: Ricard Bordalba                                               %
%                                                                         %
%   Institut de Robotica i Informatica Industrial (CSIC-UPC)              %
%   Kinematics and robot design group                                     %
%                                                                         %
%   Description: script to plot results from optimisation.                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Convergence plots                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_convergence
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.4, 0.04, 0.3, 0.3]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

    nexttile;
    hold on
    for optMode = optModeList % optmode
        optModeIndex = find(optModeList==optMode);
        plot(s.stats{optModeIndex,d_plot}.iterations.obj,'LineWidth',2)
        xlabel('Iteration number','FontSize',14)
        title('Cost','FontSize',14)
        %set(gca,'YScale','log')
    end
    ylim([1.5 12])
    legend(optModeList)
    set(gca,'FontSize',14)
    box on
    exportgraphics(gcf,'convergence.png','Resolution',300)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Kinematic plots                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_kinematic_error
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.88]);
    %set(gcf, 'OuterPosition', [0, 57,  1024, 1368]);
    a = [];
    t = tiledlayout(length(optModeList),1,'TileSpacing','Compact','Padding','Compact');
    for optMode = optModeList % optmode
        optModeIndex = find(optModeList==optMode);
        a(optModeIndex) = nexttile;
        box on
        hold on
        yline(0,'color',0.6*[1 1 1]);

        % plot knot point prime in projection method
        plot(s.ts_kinematic_error{optModeIndex,d_plot},'b','Linewidth',0.5);       
        if strcmp(optMode,"Projection")
            for r=1:length(s.knotTimePrime{optModeIndex,d_plot})-1
                ccc = plot(resample(s.ts_kinematic_error{optModeIndex,d_plot},[s.knotTimePrime{optModeIndex,d_plot}(r), s.knotTime{optModeIndex,d_plot}(r+1)]),'r','LineWidth',1.5);
            end
        end
        c = plot(resample(s.ts_kinematic_error{optModeIndex,d_plot},s.knotTime{optModeIndex,d_plot}),'ko','MarkerSize',3.5,'MarkerFaceColor','k');

        ylabel('ek','FontSize',12);%,'Interpreter','latex');
        xlim([0. t_f])
        box on
        ylim([-0.05 0.12])
        ylim([1e-18 1])
        yticks([1e-15 1e-10 1e-5 1])
        set(gca,'YScale','log')
        set(gca,'FontSize',11)
        xlabel('t [s]','FontSize',11);%,'Interpreter','latex');
        title('Kinematic error (' + optMode + ' method)','FontSize',11);
    end

    linkaxes(a)
    exportgraphics(gcf,'kinematic_error.png','Resolution',300)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Dynamic error plots                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_dynamic_error
    figure; 
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.3, 0.95]);
    t = tiledlayout(length(optModeList),1,'TileSpacing','Compact','Padding','Compact');
    for optMode = optModeList % optmode
        optModeIndex = find(optModeList==optMode);
        a(optModeIndex) = nexttile;
        hold on
        yline(0,'color',0.6*[1 1 1]);

        plot(s.ts_dynamic_error{optModeIndex,d_plot},'b','LineWidth',2);       
        if ~strcmp(optMode,"PKT")
            c = plot(resample(s.ts_dynamic_error{optModeIndex,d_plot},s.knotTime{optModeIndex,d_plot}),'ko','MarkerSize',6,'MarkerFaceColor','k');%,'k.','MarkerSize',20);   
        end
        plot(resample(s.ts_dynamic_error{optModeIndex,d_plot},s.colTime{optModeIndex,d_plot}),'ro','MarkerSize',6,'MarkerFaceColor','r');%'r.','MarkerSize',20);
        title('Dynamic error (' + optMode + ' method)');
        xlim([0 s.ts_dynamic_error{optModeIndex,d_plot}.Time(end)])
        xlabel('t [s]');
        ylabel('ed')
        xlim([3.2 3.7])
        ylim([1e-16 250])
        yticks([1e-15 1e-10 1e-5 1e0])
        box on
        set(gca,'FontSize',14)
        set(gca,'YScale','log')
    end
    linkaxes(a)
    exportgraphics(gcf,'dynamic_error.png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            State-action plots                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_state_and_action_trajectory
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.85]);
    t = tiledlayout(length(optModeList)+1,2,'TileSpacing','Compact','Padding','Compact');
    a = [];
    % reference
    nexttile;
    box on
    hold on
    % plot configuration and knot points
    plot(ts_q,'Linewidth',0.5);   
    title('Configuration profile (Init)','FontSize',16);
    xlim([0. t_f])
    xlabel('t','FontSize',12);
    box on

    nexttile
    box on
    hold on
    % plot configuration and knot points
    plot(ts_u,'b','Linewidth',0.5);       

    title('Input profile (init)','FontSize',16);
    xlim([0. t_f])
    ylim([-max(p.umax) max(p.umax)])
    xlabel('t','FontSize',12);
    box on

    for optMode = optModeList % optmode
        a(optModeIndex) = nexttile;
        box on
        hold on
        optModeIndex = find(optModeList==optMode);
        % plot configuration and knot points
        plot(s.ts_q_opt{optModeIndex,d_plot},'Linewidth',0.5);       
        plot(resample(s.ts_q_opt{optModeIndex,d_plot},s.knotTime{optModeIndex,d_plot}),'o','MarkerSize',4);%,'MarkerFaceColor','k');
        title('Configuration profile (' + optMode + ')','FontSize',16);
        xlim([0. t_f])
        xlabel('t','FontSize',12);
        box on

        nexttile
        box on
        hold on
        % plot configuration and knot points
        plot(s.ts_u_opt{optModeIndex,d_plot},'b','Linewidth',0.5);       
        plot(resample(s.ts_u_opt{optModeIndex,d_plot},s.knotTime{optModeIndex,d_plot}),'ko','MarkerSize',4);%,'MarkerFaceColor','k');

        title('Input profile (' + optMode + ')','FontSize',16);
        xlim([0. t_f])
        ylim([-max(p.umax) max(p.umax)])
        xlabel('t','FontSize',12);
        box on
    end
    linkaxes(a)
    exportgraphics(gcf,'solution_trajectory.png')
end

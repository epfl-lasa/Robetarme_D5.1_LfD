function [ output_args ] = set_figure_for_two_clum(f,plot_options)
%-----------------------------------------
% example to use this function:
% 
% [f]=plot_vel_force(Data_pos_vel_fz,vel_samples,vel_size);
% plot_options.x_name='$x_1$ [m]';
% plot_options.y_name='$x_2$ [m]';
% plot_options.z_name='$x_3$ [m]';
% plot_options.savefilename='circleDSmodela';
% set_figure_for_two_clum(f,plot_options)
%% ----------------------------------------------

set(gca, 'TickLabelInterpreter', 'latex','FontSize',9);

if isfield(plot_options,'x_name')
    x_name = plot_options.x_name;
    xlabel(x_name, 'Interpreter','LaTex','FontSize',12)
else
end

if isfield(plot_options,'y_name')
    y_name = plot_options.y_name;
    ylabel(y_name, 'Interpreter','LaTex','FontSize',12)
else
end

if isfield(plot_options,'z_name')
    z_name = plot_options.z_name;
    zlabel(z_name, 'Interpreter','LaTex','FontSize',12)
else
end

if isfield(plot_options,'lengend_name')
    % xRange = xlim;
    yRange = ylim;
    % xlim([xRange(1) xRange(2)*1.2]);
    ylim([yRange(1) yRange(2)*1.1]);

    lengend_name = plot_options.lengend_name;
    lgd = legend(lengend_name,'Interpreter','LaTex','FontSize',10, 'Location', 'north')
    pos = get(lgd, 'Position');
    newPos = [pos(1), pos(2)*1.03, pos(3), pos(4)]; 
    set(lgd, 'Position', newPos);
    lgd.Box='off';
    lgd.Orientation='horizontal';
else
end

if isfield(plot_options,'title_name')
    title_name = plot_options.title_name;
    title(title_name, 'Interpreter','LaTex','FontSize',12)
else
end

if strcmp(get(f, 'WindowStyle'), 'docked')
    set(f, 'WindowStyle', 'normal');
end

set(gcf, 'Units', 'inches', 'Position', [5, 5, 3.75, 2.75]);

% axesPos = [0.1, 0.2, 1, 0.7]; % [left, bottom, width, height]
% set(gca, 'Position', axesPos);

set(f, 'PaperPositionMode', 'auto');

if isfield(plot_options,'savefilename')
    savefilename = plot_options.savefilename;
    print(f, [savefilename '.eps'], '-depsc', '-painters', '-r2000');  
    % print(f, [savefilename '.eps'], '-depsc',  '-opengl', '-r1000');  
%     print(f, [savefilename '.png'], '-dpng', '-r3000');  
    % print(f,[savefilename '.pdf'], '-dpdf', '-painters', '-r1000');  

else
end

end


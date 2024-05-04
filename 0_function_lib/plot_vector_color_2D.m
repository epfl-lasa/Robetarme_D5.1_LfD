function [h_data, h_vel] = plot_vector_color_2D( Data_pose_vel,vector_param )
% Data_pose_vel 123 vector is 456 nothing is 789 is position

%% use parameter
vel_samples=vector_param.vel_samples;
vel_size=vector_param.vel_size;
color=vector_param.color;

if isempty(vector_param.vel_samples)
    vel_samples=1;
end
if isempty(vector_param.vel_size)
    vel_size=1;
end
if isempty(vector_param.color)
    color=[1 0 0];
end
%%     figure %3D
% subplot(212)
    h_data = plot(Data_pose_vel(3,:),Data_pose_vel(4,:),'Color',[0 0.45 0.74],'Linewidth',2); hold on;
    h_att = [];
    att=[Data_pose_vel(3,end),Data_pose_vel(4,end)];
%     h_att = scatter(att(1),att(2),150,[0 0 0],'d','Linewidth',2); hold on;
    
    % Plot Velocities of Reference Trajectories
    vel_points = Data_pose_vel([3 4 1 2],1:vel_samples:end);
    U = zeros(size(vel_points,2),1);
    V = zeros(size(vel_points,2),1);
    for i = 1:size(vel_points, 2)
        dir_    = vel_points(3:end,i);%/norm(vel_points(4:end,i))
        U(i,1)   = dir_(1);
        V(i,1)   = dir_(2);
    end
    h_vel = quiver(vel_points(1,:)',vel_points(2,:)', U, V,...
        vel_size, 'Color', color, 'LineWidth',1); hold on;
    grid on;axis equal;
    box on;
a=1;
end


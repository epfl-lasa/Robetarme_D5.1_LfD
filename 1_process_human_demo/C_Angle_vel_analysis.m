%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rui Wu 2022.04.18
%   read human demos data and find feature from frequence field
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all; close all;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set path
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_of_data='2ndSessionRobetarmeRecording';

demo_name='processed';

exp_kind='small_plant_shot_trans_to_robot_framework';

if isunix
    %-----  path for linux
    path_of_load = ['./0_human_demo/' folder_of_data '/' demo_name '/'];
    path_of_plot=['./0_figure/' exp_kind '/'];
else
    path_of_load = ['.\0_human_demo\' folder_of_data '\' demo_name '\'];
    path_of_plot=['.\0_figure\' exp_kind '\'];
end

status = mkdir(path_of_plot);   %path for save figure

%% read data



load([path_of_load exp_kind])

pos_x=[];pos_y=[];pos_z=[];
vel_x=[];vel_y=[];vel_z=[];
ang_x=[];ang_y=[];ang_z=[];
w_x=[];w_y=[];w_z=[];
quat_w=[];quat_x=[];quat_y=[];quat_z=[];


% expname='subject1';
% for exptime=1:4

% expname='subject2 linear';
% for exptime=5:8

% expname='subject2 tryC';
% for exptime=8

expname='subject3';
for exptime=9:12

% expname='subject_all';
% for exptime=1:12

exptime
    clear Pose Vel pause_time time_squence time quaternion ang AngleVel quaternion_d X_filloutline X_filter time_diff x y
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Option1: choice one to analysis
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Pose(:,1:3)=proc_data{exptime}.X(1:3,:)';
    Vel(:,1:3)=proc_data{exptime}.X_dot(1:3,:)';
    time=proc_data{exptime}.time;
    quaternion(:,2:4)=proc_data{exptime}.Q_xyzw(1:3,:)';
    quaternion(:,1)=proc_data{exptime}.Q_xyzw(4,:)';
    AngleVel(:,1:3)=proc_data{exptime}.AngleVelocity(1:3,:)';
    ang(:,1:3)=proc_data{exptime}.Angle(1:3,:)';

    target_pose=proc_data{exptime}.target_pose(1:3,:)';
    target_Qwxyz=proc_data{exptime}.target_Qwxyz(1:4,:)';

    operation_time=time(end)-time(1);

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Use PCA
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use_pca=1;
    if use_pca==1
        %--- pac for pose
            [mappedX, mapping] = pca_drtoolbox(Pose);
    
            mappedX=[zeros([length(mappedX(:,1)),1]) mappedX];
    
            [coeff, score,latent]=pca(Pose);
            PoseAfterPca=Pose*coeff;


            VelAfterPca=Vel*coeff;
            AngleVelAfterPca=AngleVel*coeff;


    end

    %% save data in one
    pos_x = [pos_x; Pose(:,1)];
    pos_y = [pos_y; Pose(:,2)];
    pos_z = [pos_z; Pose(:,3)];

    vel_x = [vel_x; Vel(:,1)];
    vel_y = [vel_y; Vel(:,2)];
    vel_z = [vel_z; Vel(:,3)];

    quat_w = [quat_w; quaternion(:,1)]; 
    quat_x = [quat_x; quaternion(:,2)];
    quat_y = [quat_y; quaternion(:,3)];
    quat_z = [quat_z; quaternion(:,4)];

   
    w_x = [w_x; AngleVel(:,1)]; 
    w_y = [w_y; AngleVel(:,2)];
    w_z = [w_z; AngleVel(:,3)];

    ang_x = [ang_x; ang(:,1)]; 
    ang_y = [ang_y; ang(:,2)]; 
    ang_z = [ang_z; ang(:,2)]; 


    fig_out=figure('Position', [100 100 800 350]);
    subplot(221)
    plot3(Vel(:,1),Vel(:,2),Vel(:,3),'b');hold on;grid on;
    plot3(VelAfterPca(:,1),VelAfterPca(:,2),VelAfterPca(:,3),'r');hold on;grid on;
    xlabel('v_x');ylabel('v_y');zlabel('v_z')
    axis equal
    subplot(222)
    plot3(AngleVel(:,1),AngleVel(:,2),AngleVel(:,3),'b');grid on;hold on;
    plot3(AngleVelAfterPca(:,1),AngleVelAfterPca(:,2),AngleVelAfterPca(:,3),'r');grid on;
    xlabel('w_x');ylabel('w_y');zlabel('w_z')
    axis equal
    sgtitle([num2str(exptime) ': ' expname ' OPtime= ' num2str(operation_time) 's'])
    %% save crete results plot
    print(fig_out,[path_of_plot num2str(exptime) '-' expname],'-depsc','-tiff');
    print(fig_out,[path_of_plot num2str(exptime) '-' expname],'-dpng');
    !echo Save figure finished!
    
    %% boxline figure to analysis angle velocity
    % subplot(223)
    % b=boxchart(Vel);
    % b.JitterOutliers = 'on';
    % b.MarkerStyle = '.';
    % boxplot(Vel,'Whisker',5.5);grid on;
    % 
    % subplot(224)
    % b=boxchart(AngleVel);
    % b.JitterOutliers = 'on';
    % b.MarkerStyle = '.';
    % % ylim([-15 15])
    % boxplot(AngleVel,'Whisker',5.5);grid on;

    %% --- swarmchart
    subplot(223)
    x = [ones(1,length(Vel)) 2*ones(1,length(Vel)) 3*ones(1,length(Vel))];
    y=[Vel(:,1)' Vel(:,2)' Vel(:,3)'];
    b=swarmchart(x,y,1);grid on;

    subplot(224)
    x = [ones(1,length(Vel)) 2*ones(1,length(Vel)) 3*ones(1,length(Vel))];
    y=[AngleVel(:,1)' AngleVel(:,2)' AngleVel(:,3)'];
    b=swarmchart(x,y,1);grid on;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analysis data together
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analysisAll=1
if analysisAll==1
    fig_out1=figure('Position', [100 100 800 350]);
    subplot(221)
    plot3(vel_x,vel_y,vel_z);grid on;
    xlabel('v_x');ylabel('v_y');zlabel('v_z')
    axis equal
    subplot(222)
    plot3(w_x,w_y,w_z);grid on;
    xlabel('w_x');ylabel('w_y');zlabel('w_z')
    axis equal

    % %--- boxchart
    % subplot(223)
    % b=boxchart(Vel);
    % b.JitterOutliers = 'on';
    % b.MarkerStyle = '.';grid on;
    % % boxplot(Vel,'Whisker',5.5);grid on;
    % 
    % subplot(224)
    % b=boxchart(AngleVel);
    % b.JitterOutliers = 'on';
    % b.MarkerStyle = '.';grid on;
    % % ylim([-15 15])
    % % boxplot(AngleVel,'Whisker',5.5);

    %--- swarmchart
    subplot(223)
    x = [ones(1,length(vel_x)) 2*ones(1,length(vel_x)) 3*ones(1,length(vel_x))];
    y=[vel_x' vel_y' vel_z'];
    b=swarmchart(x,y,1);grid on;
    ylim([-0.5 0.5])

    subplot(224)
    x = [ones(1,length(w_x)) 2*ones(1,length(w_x)) 3*ones(1,length(w_x))];
    y=[w_x' w_y' w_z'];
    b=swarmchart(x,y,1);grid on;
    ylim([-1.5 1.5])

    sgtitle(expname)
    %% save crete results plot
    print(fig_out1,[path_of_plot expname],'-depsc','-tiff');
    print(fig_out1,[path_of_plot expname],'-dpng');
    !echo Save figure finished!
end

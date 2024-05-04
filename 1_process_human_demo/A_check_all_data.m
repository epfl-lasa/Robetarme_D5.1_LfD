%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rui Wu 2022.04.18
%   read human demos data and find feature from frequence field
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; close all;

color_r=[0.85 0.33 0.1];
color_b=[0 0.45 0.74];
color_g=[0.47 0.67 0.19];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set path
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_of_data='2ndSessionRobetarmeRecording';

demo_name='processed';

if isunix
    %-----  path for linux
    path_of_load = ['./0_human_demo/' folder_of_data '/' demo_name '/'];
    path_of_plot=['./0_figure/' folder_of_data '/'];
else
    path_of_load = ['.\0_human_demo\' folder_of_data '\' demo_name '\'];
    path_of_plot=['.\0_figure\' folder_of_data '\'];
end

status = mkdir(path_of_plot);   %path for save figure

%% read data

exp_kind='small_plant_shot_trans_to_robot_framework';

load([path_of_load exp_kind])


%% one by one
time_mean_std=[];
distance_mean_std=[];
angle_mean_std=[];
distance_mean_std_all=[];
angle_mean_std_all=[];
vel_mean_std_all=[];

for i=1:size(proc_data,2)
    clear Pose Vel pause_time Fs time_squence time quaternion_wxyz target_Qwxyz target_pose
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Option1: choice one to analysis
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pose(:,2:4)=proc_data{i}.X(1:3,:)';
    Vel(:,2:4)=proc_data{i}.X_dot(1:3,:)';
    pause_time=proc_data{i}.dt;
    time=proc_data{i}.time;
    Fs=1/pause_time(1);
    % time_squence(:,1)=linspace(0,(length(Pose(:,2))-1)*pause_time,length(Pose(:,2)));


    quaternion_wxyz(:,2:4)=proc_data{i}.Q_xyzw(1:3,:)';
    quaternion_wxyz(:,1)=proc_data{i}.Q_xyzw(4,:)';
    quaternion_xyzw=proc_data{i}.Q_xyzw(1:4,:)';

    target_pose=proc_data{i}.target_pose(1:3,:)';
    target_Qwxyz=proc_data{i}.target_Qwxyz(1:4,:)';



    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot with time
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% downsample data
    downsample=0;
    if downsample
        nb_points=3000; %downsample to how many point
        number_before_downsample=length(Pose(:,2));
        [Pose_downsample] = downsample_3d_wr(Pose(:,2:4), nb_points);
    else
        Pose_downsample=Pose(:,2:4);
    end
    
    pause_time_downsample=pause_time;
    
    %% plot trajectory
    fig_name=[folder_of_data '_circleDS_with_linearDS.gif'];
    save_fig_path=[path_of_plot fig_name];

    
    p2=figure
    subplot(3,2,[1 3 5])
    plot3(Pose_downsample(:,1),Pose_downsample(:,2),Pose_downsample(:,3),...
            '-*','Color','b','MarkerSize',1);
    grid on;
    axis equal;
    X_limit=get(gca,'Xlim'); 
    Y_limit=get(gca,'Ylim'); 
    Z_limit=get(gca,'Zlim'); 
    axis_limt=[X_limit Y_limit Z_limit];
    xlabel('x');ylabel('y');zlabel('z');


    quaternion_wxyz = quaternion_wxyz / norm(quaternion_wxyz);
    target_Qwxyz = target_Qwxyz / norm(target_Qwxyz);
    theta_deg=[];
    for j=1:length(quaternion_wxyz)
        theta_deg = [theta_deg; rad2deg(2 * acos(abs(dot(quaternion_wxyz(j,:), target_Qwxyz))))];
        vel_norm=norm(Vel(j,2:4));
    end
    time_mean_std=[time_mean_std; time(end)-time(1) ;];
    distance_mean_std=[distance_mean_std; abs(mean(Pose_downsample(:,1))-target_pose(1)) std(Pose_downsample(:,1)-target_pose(1)-(mean(Pose_downsample(:,1))-target_pose(1)));];
    angle_mean_std=[angle_mean_std; mean(theta_deg) std(theta_deg);];
    

    distance_mean_std_all{i}=[Pose_downsample(:,1)-target_pose(1)];
    angle_mean_std_all{i}=[theta_deg];
    vel_mean_std_all{i}=[vel_norm];

    title({['operation time = ' num2str(time(end)-time(1)) ' s']; ...
        ['average dis = ' num2str(abs(mean(Pose_downsample(:,1))-target_pose(1))) '±' num2str(std(Pose_downsample(:,1)-target_pose(1)-(mean(Pose_downsample(:,1))-target_pose(1)))) 'm']; ...
        ['average angle = ' num2str(mean(theta_deg)) '±' num2str(std(theta_deg)) 'degree']})
    sgtitle(['exp ' num2str(i) ' all data'])
    

    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% frequence analysis
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_data=[Pose(:,2:4)];
    
    frequence_limt=[0,5];
    frequence_limt_log=[0,10^2];
    %--- plot ori signal and do Fourier transform
    for i=1:length(plot_data(1,:))
        x=plot_data(:,i);
        
        %--- Fourier transform
        y = fft(x); 
        f = (0:length(y)-1)*Fs/length(y);%frequence
        
        %--- turn the negative
        n = length(x);                         
        fshift = (-n/2:n/2-1)*(Fs/n);
        yshift = fftshift(y);
        figure(p2)
        subplot(length(plot_data(1,:)),2,i*2)
        plot(fshift,abs(yshift));
        
        %--- different axis
        case_axis=1; 
        switch case_axis
            case 0
                plot(fshift,abs(yshift))
                xlim(frequence_limt_log)
                xlabel('Frequence');
                ylabel('Magnitude');
            case 1
                loglog(fshift,abs(yshift))
                xlim(frequence_limt_log)
                xlabel('log(Frequence)');
                ylabel('log(Magnitude)');
            case 2
                semilogx(fshift,abs(yshift))
                xlim(frequence_limt_log)
                xlabel('log(Frequence)');
                ylabel('Magnitude');
            case 3
                semilogy(fshift,abs(yshift))
                xlim(frequence_limt_log)
                xlabel('Frequence');
                ylabel('log(Magnitude)');
        end
    end
    
    
    %% plot trajectory in gif
    plot_with_time=0;
    if plot_with_time==1
        figure
        p1=animatedline(Pose_downsample(1,1),Pose_downsample(1,2),Pose_downsample(1,3));
        p1.Color = 'r';
        p1.LineWidth = 1.0;
        p1.LineStyle = '-';
        xlabel('x');ylabel('y');zlabel('z');
        legend('trj');grid on;
        title('trajectory');
        axis equal;
        axis(axis_limt);
        view(77,11);
        grid on;
    
        for i=1:length(Pose_downsample)
            %%% plot trj    
            addpoints(p1,Pose_downsample(i,1),Pose_downsample(i,2),Pose_downsample(i,3));  
            drawnow
            pause(pause_time_downsample);
            hold on;
    
        end
    end


end


%% calculate distance angle mean std for all
distance_mean_std_all_1=[];
angle_mean_std_all_1=[];
distance_mean_std_all_2=[];
angle_mean_std_all_2=[];
distance_mean_std_all_3=[];
angle_mean_std_all_3=[];
vel_mean_std_1=[];
vel_mean_std_2=[];
vel_mean_std_3=[];

for i=1:4
    distance_mean_std_all_1=[distance_mean_std_all_1;distance_mean_std_all{i}];
    angle_mean_std_all_1=[angle_mean_std_all_1;angle_mean_std_all{i}];
    vel_mean_std_1=[vel_mean_std_1;vel_mean_std_all{i}];
end
for i=5:8
    distance_mean_std_all_2=[distance_mean_std_all_2;distance_mean_std_all{i}];
    angle_mean_std_all_2=[angle_mean_std_all_2;angle_mean_std_all{i}];
    vel_mean_std_2=[vel_mean_std_2;vel_mean_std_all{i}];
end
for i=9:12
    distance_mean_std_all_3=[distance_mean_std_all_3;distance_mean_std_all{i}];
    angle_mean_std_all_3=[angle_mean_std_all_3;angle_mean_std_all{i}];
    vel_mean_std_3=[vel_mean_std_3;vel_mean_std_all{i}];
end

all_dis_angel_subject=[mean(distance_mean_std_all_1) std(distance_mean_std_all_1);
                       mean(distance_mean_std_all_2) std(distance_mean_std_all_2);
                       mean(distance_mean_std_all_3) std(distance_mean_std_all_3);]

all_dis_angel_subject=[mean(angle_mean_std_all_1) std(angle_mean_std_all_1);
                       mean(angle_mean_std_all_2) std(angle_mean_std_all_2);
                       mean(angle_mean_std_all_3) std(angle_mean_std_all_3);]

all_ve_subject=[mean(vel_mean_std_1) std(vel_mean_std_1);
                       mean(vel_mean_std_2) std(vel_mean_std_2);
                       mean(vel_mean_std_3) std(vel_mean_std_3);]
    
%% plot 
% Data from the table
subjects = {'Subject 1', 'Subject 3', 'Subject 2'};
nozzle_motion_velocity = [0.31, 0.33, 0.19];
nozzle_motion_velocity_std = [0.10, 0.20, 0.10];
nozzle_angle = 90-[89.56, 89.55, 89.76];
nozzle_angle_std = [0.15, 0.13, 0.09];
nozzle_distance = [-0.85, -0.84, -0.62];
nozzle_distance_std = [0.07, 0.06, 0.06];
rebound_rate = [7.5, 7.5, 43.6]; % NaN for missing data

% Create figure and subplots
figure;

% Nozzle motion velocity
subplot(4,1,1);
bar(nozzle_motion_velocity, 'FaceColor', [0.2, 0.2, 0.5]);
hold on;
errorbar(1:3, nozzle_motion_velocity, nozzle_motion_velocity_std, '.k');
title('Nozzle Motion Velocity [m/s]');
set(gca, 'xticklabel', subjects, 'XTick',1:numel(subjects));
ylabel('Velocity [m/s]');

% Nozzle angle
subplot(4,1,2);
bar(nozzle_angle, 'FaceColor', [0.2, 0.5, 0.2]);
hold on;
errorbar(1:3, nozzle_angle, nozzle_angle_std, '.k');
title('Nozzle Angle [degree]');
set(gca, 'xticklabel', subjects, 'XTick',1:numel(subjects));
ylabel('Angle [degree]');

% Nozzle distance
subplot(4,1,3);
bar(nozzle_distance, 'FaceColor', [0.5, 0.2, 0.2]);
hold on;
errorbar(1:3, nozzle_distance, nozzle_distance_std, '.k');
title('Nozzle Distance [m]');
set(gca, 'xticklabel', subjects, 'XTick',1:numel(subjects));
ylabel('Distance [m]');

% Rebound rate
subplot(4,1,4);
bar(rebound_rate, 'FaceColor', [0.5, 0.5, 0.2]); % Only Subject 2 has a value
hold on;
% plot([1, 1], rebound_rate_range, '-r', 'LineWidth', 2); % Range for Subject 1
title('Rebound Rate');
set(gca, 'xticklabel', subjects, 'XTick',1:numel(subjects));
ylabel('Rate [%]');

% Adjust layout
set(gcf, 'Position', [100, 100, 600, 800]); % Resize figure

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot all with target
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos_x=[];pos_y=[];pos_z=[];
vel_x=[];vel_y=[];vel_z=[];
ang_x=[];ang_y=[];ang_z=[];
w_x=[];w_y=[];w_z=[];
quat_w=[];quat_x=[];quat_y=[];quat_z=[];



for exptime=1:size(proc_data,2)

% expname='subject2 linear';
% for exptime=5:7

% expname='subject2 tryC';
% for exptime=8

% expname='subject3';
% for exptime=9:12

exptime
    clear Pose Vel pause_time time_squence time quaternion_wxyz ang AngleVel quaternion_d X_filloutline X_filter time_diff
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Option1: choice one to analysis
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Pose(:,1:3)=proc_data{exptime}.X(1:3,:)';
    Vel(:,1:3)=proc_data{exptime}.X_dot(1:3,:)';
    time=proc_data{exptime}.time;
    quaternion_wxyz(:,2:4)=proc_data{exptime}.Q_xyzw(1:3,:)';
    quaternion_wxyz(:,1)=proc_data{exptime}.Q_xyzw(4,:)';
    AngleVel(:,1:3)=proc_data{exptime}.AngleVelocity(1:3,:)';
    ang(:,1:3)=proc_data{exptime}.Angle(1:3,:)';

    target_pose=proc_data{exptime}.target_pose(1:3,:)';
    target_Qwxyz=proc_data{exptime}.target_Qwxyz(1:4,:)';

    operation_time=time(end)-time(1);

    pos_x = [pos_x; Pose(:,1)];
    pos_y = [pos_y; Pose(:,2)];
    pos_z = [pos_z; Pose(:,3)];

    vel_x = [vel_x; Vel(:,1)];
    vel_y = [vel_y; Vel(:,2)];
    vel_z = [vel_z; Vel(:,3)];

    quat_w = [quat_w; quaternion_wxyz(:,1)]; 
    quat_x = [quat_x; quaternion_wxyz(:,2)];
    quat_y = [quat_y; quaternion_wxyz(:,3)];
    quat_z = [quat_z; quaternion_wxyz(:,4)];

   
    w_x = [w_x; AngleVel(:,1)]; 
    w_y = [w_y; AngleVel(:,2)];
    w_z = [w_z; AngleVel(:,3)];

    ang_x = [ang_x; ang(:,1)]; 
    ang_y = [ang_y; ang(:,2)]; 
    ang_z = [ang_z; ang(:,2)]; 


    fig = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
    for i=1:3
        subplot(4,5,1+(i-1)*5)
        plot(time,ang(:,i))
        title(['ang' num2str(i)])
    end
    for i=1:3
        subplot(4,5,2+(i-1)*5)
        plot(time,AngleVel(:,i))
        title(['AngleVel' num2str(i)])
    end
    
    for i=1:3
        subplot(4,5,3+(i-1)*5)
        plot(time,Pose(:,i))
        title(['Pose' num2str(i)])
    end
    for i=1:4
        subplot(4,5,5+(i-1)*5)
        plot(time,quaternion_wxyz(:,i))
        title(['quaternion' num2str(i)])
    end
    for i=1:3
        subplot(4,5,4+(i-1)*5)
        plot(time,Vel(:,i))
        title(['Vel' num2str(i)])
    end

    %--- plot trj
    subplot(4,4,[13:15])
    plot3(Pose(:,1),Pose(:,2),Pose(:,3));hold on;grid on;
    xlabel('x');ylabel('y');zlabel('z');axis equal
    % xlim([-2 2]);ylim([-2 2]);zlim([1 2])
    %--- plot target surface
    point1=target_pose+[ 0 0.25 0.25];point2=target_pose+[ 0 -0.25 0.25];
    point3=target_pose+[ 0 -0.25 -0.25];point4=target_pose+[ 0 0.25 -0.25];

    X_patch=[point1(1) point2(1) point3(1) point4(1)];
    Y_patch=[point1(2) point2(2) point3(2) point4(2)];
    Z_patch=[point1(3) point2(3) point3(3) point4(3)];
    patch(X_patch,Y_patch,Z_patch,[0 .5 .5])

    sgtitle(['exp ' num2str(exptime) ' all data'])


end


for exptime=1:size(proc_data,2)

% expname='subject2 linear';
% for exptime=5:7

% expname='subject2 tryC';
% for exptime=8

% expname='subject3';
% for exptime=9:12

exptime
    clear Pose Vel pause_time time_squence time quaternion_wxyz ang AngleVel quaternion_d X_filloutline X_filter time_diff
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Option1: choice one to analysis
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Pose(:,1:3)=proc_data{exptime}.X(1:3,:)';
    Vel(:,1:3)=proc_data{exptime}.X_dot(1:3,:)';
    time=proc_data{exptime}.time;
    quaternion_wxyz(:,2:4)=proc_data{exptime}.Q_xyzw(1:3,:)';
    quaternion_wxyz(:,1)=proc_data{exptime}.Q_xyzw(4,:)';
    AngleVel(:,1:3)=proc_data{exptime}.AngleVelocity(1:3,:)';
    ang(:,1:3)=proc_data{exptime}.Angle(1:3,:)';

    target_pose=proc_data{exptime}.target_pose(1:3,:)';
    target_Qwxyz=proc_data{exptime}.target_Qwxyz(1:4,:)';

    operation_time=time(end)-time(1);

    pos_x = [pos_x; Pose(:,1)];
    pos_y = [pos_y; Pose(:,2)];
    pos_z = [pos_z; Pose(:,3)];

    vel_x = [vel_x; Vel(:,1)];
    vel_y = [vel_y; Vel(:,2)];
    vel_z = [vel_z; Vel(:,3)];

    quat_w = [quat_w; quaternion_wxyz(:,1)]; 
    quat_x = [quat_x; quaternion_wxyz(:,2)];
    quat_y = [quat_y; quaternion_wxyz(:,3)];
    quat_z = [quat_z; quaternion_wxyz(:,4)];

   
    w_x = [w_x; AngleVel(:,1)]; 
    w_y = [w_y; AngleVel(:,2)];
    w_z = [w_z; AngleVel(:,3)];

    ang_x = [ang_x; ang(:,1)]; 
    ang_y = [ang_y; ang(:,2)]; 
    ang_z = [ang_z; ang(:,2)]; 


    fig = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
    for i=1:3
        subplot(4,5,1+(i-1)*5)
        plot(time,ang(:,i))
        title(['ang' num2str(i)])
    end
    for i=1:3
        subplot(4,5,2+(i-1)*5)
        plot(time,AngleVel(:,i))
        title(['AngleVel' num2str(i)])
    end
    
    for i=1:3
        subplot(4,5,3+(i-1)*5)
        plot(time,Pose(:,i))
        title(['Pose' num2str(i)])
    end
    for i=1:4
        subplot(4,5,5+(i-1)*5)
        plot(time,quaternion_wxyz(:,i))
        title(['quaternion' num2str(i)])
    end
    for i=1:3
        subplot(4,5,4+(i-1)*5)
        plot(time,Vel(:,i))
        title(['Vel' num2str(i)])
    end

    %--- plot trj
    subplot(4,4,[13:15])
    plot3(Pose(:,1),Pose(:,2),Pose(:,3));hold on;grid on;
    xlabel('x');ylabel('y');zlabel('z');axis equal
    %--- plot target surface
    point1=target_pose+[ 0 0.25 0.25];point2=target_pose+[ 0 -0.25 0.25];
    point3=target_pose+[ 0 -0.25 -0.25];point4=target_pose+[ 0 0.25 -0.25];

    X_patch=[point1(1) point2(1) point3(1) point4(1)];
    Y_patch=[point1(2) point2(2) point3(2) point4(2)];
    Z_patch=[point1(3) point2(3) point3(3) point4(3)];
    patch(X_patch,Y_patch,Z_patch,[0 .5 .5])
    sgtitle(['exp ' num2str(exptime) ' all data after change'])

end

%% save mat file
save_after_change=0;
if save_after_change==1
save([path_of_load exp_kind],'proc_data');
!echo SAVE TO MAT DONE!
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Option2: merges all to analysis
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% time analysis
for i=1:size(proc_data,2)
    time_all(i)=proc_data{i}.time(end)-proc_data{i}.time(1);
end
time_all'

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot for paper
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
for exptime=[7 10]
    exptime
    clear Pose Vel pause_time time_squence time quaternion_wxyz ang AngleVel quaternion_d X_filloutline X_filter time_diff
    
    Pose(:,1:3)=proc_data{exptime}.X(1:3,:)';
    Vel(:,1:3)=proc_data{exptime}.X_dot(1:3,:)';
    time=proc_data{exptime}.time;
    quaternion_wxyz(:,2:4)=proc_data{exptime}.Q_xyzw(1:3,:)';
    quaternion_wxyz(:,1)=proc_data{exptime}.Q_xyzw(4,:)';
    AngleVel(:,1:3)=proc_data{exptime}.AngleVelocity(1:3,:)';
    ang(:,1:3)=proc_data{exptime}.Angle(1:3,:)';

    target_pose=proc_data{exptime}.target_pose(1:3,:)';
    target_Qwxyz=proc_data{exptime}.target_Qwxyz(1:4,:)';

    operation_time=time(end)-time(1);

    if exptime==7
        subplot(121)
        data_num=[2000:4000]
        data_num=[1:5000]
        plot3(Pose(data_num,1),Pose(data_num,2),Pose(data_num,3),'-*','Color','b','MarkerSize',1);hold on;grid on;
    else
        data_num=[2000:3000]
        data_num=[1:5000]
        subplot(122)
        plot3(Pose(data_num,1),Pose(data_num,2),Pose(data_num,3),'-*','Color','b','MarkerSize',1);hold on;grid on;
    end
    xlabel('x');ylabel('y');zlabel('z');axis equal
end


%% for paper less data
for exptime=[7 10]
    exptime
    clear Pose Vel pause_time time_squence time quaternion_wxyz ang AngleVel quaternion_d X_filloutline X_filter time_diff
    
    Pose(:,1:3)=proc_data{exptime}.X(1:3,:)';
    Vel(:,1:3)=proc_data{exptime}.X_dot(1:3,:)';
    time=proc_data{exptime}.time;
    quaternion_wxyz(:,2:4)=proc_data{exptime}.Q_xyzw(1:3,:)';
    quaternion_wxyz(:,1)=proc_data{exptime}.Q_xyzw(4,:)';
    AngleVel(:,1:3)=proc_data{exptime}.AngleVelocity(1:3,:)';
    ang(:,1:3)=proc_data{exptime}.Angle(1:3,:)';

    target_pose=proc_data{exptime}.target_pose(1:3,:)';
    target_Qwxyz=proc_data{exptime}.target_Qwxyz(1:4,:)';

    operation_time=time(end)-time(1);

    if exptime==7
        fig1 = figure
        data_num=[2000:4000]
        plot3(Pose(data_num,1),Pose(data_num,2),Pose(data_num,3),'-*','Color',color_b,'MarkerSize',1);hold on;grid on;
        axis equal
        Fig4_a_1=Pose(data_num,1:3);
        view(-60,17)
    else
        data_num=[2000:3000]
        fig2 = figure
        plot3(Pose(data_num,1),Pose(data_num,2),Pose(data_num,3),'-*','Color',color_b,'MarkerSize',1);hold on;grid on;
        axis equal
        Fig4_a_2=Pose(data_num,1:3);
        view(-60,17)
    end
    xlabel('x');ylabel('y');zlabel('z');axis equal
end



plot_options.x_name='$x_1$ [m]';
plot_options.y_name='$x_2$ [m]';
plot_options.z_name='$x_3$ [m]';
plot_options.savefilename='demo_trj1';

set_figure_for_two_clum(fig1,plot_options)

plot_options.x_name='$x_1$ [m]';
plot_options.y_name='$x_2$ [m]';
plot_options.z_name='$x_3$ [m]';
plot_options.savefilename='demo_trj2';

set_figure_for_two_clum(fig2,plot_options)



%% for paper less data 2D
for exptime=[7 10]
    exptime
    clear Pose Vel pause_time time_squence time quaternion_wxyz ang AngleVel quaternion_d X_filloutline X_filter time_diff
    
    Pose(:,1:3)=proc_data{exptime}.X(1:3,:)';
    Vel(:,1:3)=proc_data{exptime}.X_dot(1:3,:)';
    time=proc_data{exptime}.time;
    quaternion_wxyz(:,2:4)=proc_data{exptime}.Q_xyzw(1:3,:)';
    quaternion_wxyz(:,1)=proc_data{exptime}.Q_xyzw(4,:)';
    AngleVel(:,1:3)=proc_data{exptime}.AngleVelocity(1:3,:)';
    ang(:,1:3)=proc_data{exptime}.Angle(1:3,:)';

    target_pose=proc_data{exptime}.target_pose(1:3,:)';
    target_Qwxyz=proc_data{exptime}.target_Qwxyz(1:4,:)';

    operation_time=time(end)-time(1);

    if exptime==7
        fig1 = figure
        data_num=[2000:4000]
        plot(Pose(data_num,2),Pose(data_num,3),'-*','Color',color_b,'MarkerSize',1);hold on;grid on;
        axis equal
        Fig4_a_1=Pose(data_num,1:3);
        plot_options.x_name='$y$ [m]';
        plot_options.y_name='$z$ [m]';
        plot_options.savefilename='demo_trj1';
        
        set_figure_for_two_clum(fig1,plot_options)
    else
        data_num=[2000:3000]
        fig2 = figure
        plot(Pose(data_num,2),Pose(data_num,3),'-*','Color',color_b,'MarkerSize',1);hold on;grid on;
        axis equal
        Fig4_a_2=Pose(data_num,1:3);
        plot_options.x_name='$y$ [m]';
        plot_options.y_name='$z$ [m]';
        plot_options.savefilename='demo_trj2';
        
        set_figure_for_two_clum(fig2,plot_options)
    end
end









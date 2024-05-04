%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rui Wu 2022.07.24
%   read human demos data and segment into peridico and nonpredico part
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all; close all;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set path
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_of_data='2ndSessionRobetarmeRecording';

demo_name='processed';

if isunix
    %-----  path for linux
    path_of_load = ['./0_human_demo/' folder_of_data '/processed/'];
    path_of_plot=['./0_figure/' folder_of_data '/'];
    path_of_save = ['./0_human_demo/' folder_of_data '/processed/'];
else
    path_of_load = ['.\0_human_demo\' folder_of_data '\processed\'];
    path_of_plot=['.\0_figure\' folder_of_data '\'];
    path_of_save = ['.\0_human_demo\' folder_of_data '\processed\'];
end

status = mkdir(path_of_plot);   %path for save figure

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read data
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_kind='small_plant_shot_trans_to_robot_framework_seprated';

load([path_of_load exp_kind])

subject1_p=[];subject1_np=[];subject1_time=[];
subject2_p=[];subject2_np=[];subject2_time=[];
subject3_p=[];subject3_np=[];subject3_time=[];

choiceDataSetMode=3;

if choiceDataSetMode==1
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% for learnCycleDSFromHumanDemo opt1: Xstored
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for trial_num=1:size(proc_data_seprated,2)
    
        clear data_p data_np time
    
        data_p(:,:)=proc_data_seprated{trial_num}.periodic(:,:);
        data_np(:,:)=proc_data_seprated{trial_num}.non_periodic(:,:);
        time(:,:)=proc_data_seprated{trial_num}.time(:,:);
    
        %--- calculate velocity
        for i=1:3
            data_p(:,i+3)=gradient(data_p(:,i),time);
        end
        for i=1:3
            data_np(:,i+3)=gradient(data_np(:,i),time);
        end

        
    
        if 1<=trial_num&trial_num<=4
            subject1_p=[subject1_p;data_p(:,:)];
            subject1_np=[subject1_np;data_np(:,:)];
        elseif 5<=trial_num&trial_num<=8
            subject2_p=[subject2_p;data_p(:,:)];
            subject2_np=[subject2_np;data_np(:,:)];
        elseif 9<=trial_num&trial_num<=9
            subject3_p=[subject3_p;data_p(:,:)];
            subject3_np=[subject3_np;data_np(:,:)];
        end
    
    end
    
    
    
    %% choice subject to analysis
    figure
    plot3(subject1_p(:,1),subject1_p(:,2),subject1_p(:,3),'r');hold on;grid on;
    % plot3(subject2_p(:,1)+0.02,subject2_p(:,2),subject2_p(:,3),'b');hold on;
    plot3(subject3_p(:,1)+0.04,subject3_p(:,2),subject3_p(:,3),'g');hold on;
    legend('subject1','subject2','subject3')
    
    %% save data set for trainning
    %--- downsample
%     subject3_p_downsample(:,1:3)=downsample_3d_wr(subject3_p(:,1:3),3000);
%     subject3_p_downsample(:,4:6)=downsample_3d_wr(subject3_p(:,4:6),3000);

    subject3_p_downsample(:,1:3)=subject3_p(4001:4100,1:3);
    subject3_p_downsample(:,4:6)=subject3_p(4001:4100,4:6);

    subject3_p_turn(:,1)=subject3_p_downsample(:,3);
    subject3_p_turn(:,2)=subject3_p_downsample(:,2);
    subject3_p_turn(:,3)=subject3_p_downsample(:,1);
    
    subject3_p_turn(:,4)=subject3_p_downsample(:,6);
    subject3_p_turn(:,5)=subject3_p_downsample(:,5);
    subject3_p_turn(:,6)=subject3_p_downsample(:,4);

    figure
    plot_vel(subject3_p_downsample',1,1);
    title('subject3 downsample')


    Xstored=subject3_p_downsample(:,:);
    demo_kind='subject3_all_perdico_data';
    save([path_of_save demo_kind],'Xstored');
    !echo Save data done!!

elseif choiceDataSetMode==2
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% for learnCycleDSFromHumanDemo opt1: traj
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for trial_num=1:size(proc_data_seprated,2)
    % for trial_num=10
    
        clear data_p data_np time
    
        data_p(:,:)=proc_data_seprated{trial_num}.periodic(:,:);
        data_np(:,:)=proc_data_seprated{trial_num}.non_periodic(:,:);
        time(:,:)=proc_data_seprated{trial_num}.time(:,:)';
    
        if 1<=trial_num&trial_num<=4
            subject1_p=[subject1_p;data_p(:,:)];
            subject1_np=[subject1_np;data_np(:,:)];
            subject1_time=[subject1_time;time(:,:)];
        elseif 5<=trial_num&trial_num<=8
            subject2_p=[subject2_p;data_p(:,:)];
            subject2_np=[subject2_np;data_np(:,:)];
            subject2_time=[subject2_time;time(:,:)];
        elseif 9<=trial_num&trial_num<=9
            subject3_p=[subject3_p;data_p(:,:)];
            subject3_np=[subject3_np;data_np(:,:)];
            subject3_time=[subject3_time;time(:,:)];
        end
    
    end

    %% choice subject to analysis
    figure
    plot3(subject1_p(:,1),subject1_p(:,2),subject1_p(:,3),'r');hold on;grid on;
    plot3(subject2_p(:,1)+0.02,subject2_p(:,2),subject2_p(:,3),'b');hold on;
    plot3(subject3_p(:,1)+0.04,subject3_p(:,2),subject3_p(:,3),'g');hold on;
    legend('subject1','subject2','subject3')
    
    %% save data set for trainning   
    subject3_p_downsample(:,1:3)=subject3_p(2001:2100,1:3);
    subject3_time_downsample(:,1)=subject3_time(2001:2100,1);
    
    figure
    plot3(subject3_p_downsample(:,1),subject3_p_downsample(:,2),subject3_p_downsample(:,3),'r');hold on;grid on;
    title('subject3 downsample')
    
    trial.position_real.x=subject3_p_downsample(:,1);
    trial.position_real.y=subject3_p_downsample(:,2);
    trial.position_real.z=subject3_p_downsample(:,3);
    trial.position_real.time=subject3_time_downsample(:,1);
    demo_kind='subject3_all_perdico_data';
    save([path_of_save demo_kind],'trial');
    !echo Save data done!!
    
elseif choiceDataSetMode==3
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% for learnCycleDSFromHumanDemo opt1: traj
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for trial_num=1:size(proc_data_seprated,2)
    
        clear data_p data_np time
    
        data_p(:,:)=proc_data_seprated{trial_num}.periodic(:,:);
        data_np(:,:)=proc_data_seprated{trial_num}.non_periodic(:,:);
        time(:,:)=proc_data_seprated{trial_num}.time(:,:)'-proc_data_seprated{trial_num}.time(1,1);

        for i=1:3
            data_np(:,i+3)=gradient(data_np(:,i),time);
        end
    
        if 1<=trial_num&trial_num<=4
            subject1_p=[subject1_p;data_p(:,:)];
            subject1_np=[subject1_np;data_np(:,:)];
            subject1_time=[subject1_time;time(:,:)];
        elseif 8<=trial_num&trial_num<=8
            subject2_p=[subject2_p;data_p(:,:)];
            subject2_np=[subject2_np;data_np(:,:)];
            subject2_time=[subject2_time;time(:,:)];
        elseif 9<=trial_num&trial_num<=12
            subject3_p=[subject3_p;data_p(:,:)];
            subject3_np=[subject3_np;data_np(:,:)];
            subject3_time=[subject3_time;time(:,:)];
        end
    
    end

    %% choice subject to analysis
    figure
    plot3(subject1_p(:,1),subject1_p(:,2),subject1_p(:,3),'r');hold on;grid on;
    % plot3(subject2_p(:,1)+0.02,subject2_p(:,2),subject2_p(:,3),'b');hold on;
    plot3(subject3_p(:,1)+0.04,subject3_p(:,2),subject3_p(:,3),'g');hold on;
    legend('subject1','subject2','subject3')
    
    %% save subject1 data set for trainning   
    subject1_p_downsample(:,1:3)=subject1_p(2051:3050,1:3);
    subject1_time_downsample(:,1)=subject1_time(2051:3050,1);
    
    figure
    plot3(subject1_p_downsample(:,1),subject1_p_downsample(:,2),subject1_p_downsample(:,3),'r');hold on;grid on;axis equal;
    title('subject1 downsample')
    
    traj.data=subject1_p_downsample(:,1:3);
    traj.time=subject1_time_downsample(:,1);
    traj.size=size(subject1_p_downsample(:,1:3));
    
    demo_kind='subject1_all_perdico_data';
    save([path_of_save demo_kind],'traj');
    !echo Save data done!!

    %% save subject2 data set for trainning   
    % subject2_p_downsample(:,1:3)=subject2_p(2051:3050,1:3);
    % subject2_time_downsample(:,1)=subject2_time(2051:3050,1);
    % 
    % figure
    % plot3(subject2_p_downsample(:,1),subject2_p_downsample(:,2),subject2_p_downsample(:,3),'r');hold on;grid on;axis equal;
    % title('subject2 downsample')
    % 
    % traj.data=subject2_p_downsample(:,1:3);
    % traj.time=subject2_time_downsample(:,1);
    % traj.size=size(subject2_p_downsample(:,1:3));
    % 
    % demo_kind='subject2_all_perdico_data';
    % save([path_of_save demo_kind],'traj');
    % !echo Save data done!!

    %% save subject3 perdic data set for trainning   
    subject3_p_downsample(:,1:3)=subject3_p(2051:3050,1:3);
    subject3_time_downsample(:,1)=subject3_time(2051:3050,1);
    
    figure
    plot3(subject3_p_downsample(:,1),subject3_p_downsample(:,2),subject3_p_downsample(:,3),'r');hold on;grid on;axis equal;
    title('subject3 downsample')
    
    traj.data=subject3_p_downsample(:,1:3);
    traj.time=subject3_time_downsample(:,1);
    traj.size=size(subject3_p_downsample(:,1:3));
    
    demo_kind='subject3_all_perdico_data';
    save([path_of_save demo_kind],'traj');
    !echo Save data done!!

    %% save subject3 non perdic data set for trainning   
    subject3_np_downsample=[];subject3_time_downsample=[];
    subject3_np_downsample(:,1:3)=subject3_np(1051:3050,1:3);
    subject3_time_downsample(:,1)=subject3_time(1051:3050,1);
    
    figure
    plot3(subject3_np_downsample(:,1),subject3_np_downsample(:,2),subject3_np_downsample(:,3),'r');hold on;grid on;axis equal;
    title('subject3 downsample non_perdico')
    
    traj.data=subject3_np_downsample(:,1:3);
    traj.time=subject3_time_downsample(:,1);
    traj.size=size(subject3_np_downsample(:,1:3));
    
    demo_kind='subject3_all_non_perdico_data';
    save([path_of_save demo_kind],'traj');
    !echo Save data done!!

end

































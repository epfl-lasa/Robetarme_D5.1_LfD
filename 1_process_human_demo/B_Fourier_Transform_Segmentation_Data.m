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
    path_of_load = ['./0_human_demo/' folder_of_data '/' demo_name '/'];
    path_of_plot=['./0_figure/' folder_of_data '/'];
    path_of_save = ['./0_human_demo/' folder_of_data '/processed/'];
else
    path_of_load = ['.\0_human_demo\' folder_of_data '\' demo_name '\'];
    path_of_plot=['.\0_figure\' folder_of_data '\'];
    path_of_save = ['.\0_human_demo\' folder_of_data '\processed\'];
end

status = mkdir(path_of_plot);   %path for save figure

%% read data

exp_kind='small_plant_shot_trans_to_robot_framework';

load([path_of_load exp_kind])

use_pca=0;
plot_all=0;
use_fooof_get_centerAndBw=1;
do_forier=1;
plot_segment=1;

save_or_not=0;


% for trial_num=1:size(proc_data,2)
for trial_num=1

    if trial_num==[5:8]
        use_fooof_get_centerAndBw=0;
    end

    clear Pose Vel pause_time Fs time_squence time AngleVel
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  choice one to analysis
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pose(:,2:4)=proc_data{trial_num}.X(1:3,:)';
    Vel(:,2:4)=proc_data{trial_num}.X_dot(1:3,:)';
    AngleVel(:,1:3)=proc_data{trial_num}.AngleVelocity(1:3,:)';

    pause_time=proc_data{trial_num}.dt;
    time=proc_data{trial_num}.time;
    Fs=1/mean(pause_time);
    

     
    plot_data=[Pose(:,2:4)];


    frequence_limt=[0,5];
    frequence_limt_log=[0,10];

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Use PCA
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if use_pca==1
        [mappedX, mapping] = pca_drtoolbox(plot_data);

        mappedX=[zeros([length(mappedX(:,1)),1]) mappedX];

        [coeff, score,latent]=pca(plot_data);
        dataAfterPca=plot_data*coeff;

        plot_data=[ dataAfterPca(:,3) dataAfterPca(:,1) dataAfterPca(:,2)]
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot with time
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% downsample data

    if plot_all==1
        downsample=0;
        if downsample
            nb_points=3000; %downsample to how many point
            number_before_downsample=length(Pose(:,2));
            [Pose_downsample] = downsample_3d_wr(Pose(:,2:4), nb_points);
        else
            Pose_downsample=Pose(:,2:4);
        end
        
        pause_time_downsample=pause_time;
    
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% fooof, based on the nature paper, to find best period frequency
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if use_fooof_get_centerAndBw==1
        h=figure
        for i=1:3
            dt=pause_time;% sampling interval
            N=length(plot_data(:,1));% sampling points
            t=time';% sampling time
            Fs=Fs;%Sampling frequency, which is the reciprocal of the sampling interval
            n=0:1:N-1;
            f=(Fs/N).*n;%The frequency corresponding to each point on the X axis
            x=plot_data(:,i); %signal
    
            y=fft(x);% Fourier transform to get a complex number
            
            Ay=abs(y);% modulo
            Ayy=Ay*2/N;% converted to actual amplitude
    
            
            subplot(3,3,0*3+i)
            plot(t,x);hold on;grid on;
            if i==1
                title('X direction')
            elseif i==2
                title('Y direction')
            elseif i==3
                title('Z direction')
            end
            
            
            
            freqs           =  f         %row vector of frequency values
            power_spectrum  =  abs(y')         %row vector of power values
            f_range         =  [0, 10]         %fitting range (Hz)
            %--- settings :    %fooof model settings, in a struct, including:
            settings.peak_width_limts= [2.0 3.0];
            settings.max_n_peaks= 1;
            settings.min_peak_height=0.2;
            settings.peak_threshold=2.0;
            settings.aperiodic_mode='fixed';
            return_model    = 1         %boolean of whether to return the FOOOF model fit, optional
            
            fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model)
    
            subplot(3,3,1*3+i)
            fooof_plot_subplot(h,fooof_results);hold on;grid on
            subplot(3,3,2*3+i)
            plot(fooof_results.freqs,(fooof_results.fooofed_spectrum-fooof_results.ap_fit));hold on;grid on
    
            %--- invers fuliye
            di=ifft(fooof_results.ap_fit);
            yi=ifft(fooof_results.fooofed_spectrum-fooof_results.ap_fit);
    
            plot_data_p_fooof(:,i)=real(yi);
            plot_data_np_fooof(:,i)=real(di);

            %--- manually set some peak
            if (trial_num==1)
                if i==1
                    fooof_results.peak_params(1)=2.2;
                    fooof_results.peak_params(3)=0.7;
                end
            elseif (trial_num==2)
                if i==1
                    fooof_results.peak_params(1)=2.3;
                    fooof_results.peak_params(3)=0.5;
                end
            elseif (trial_num==4)
                if i==1
                    fooof_results.peak_params(1)=2.7;
                    fooof_results.peak_params(2)=5.4;
                    fooof_results.peak_params(3)=0.4;
                end
            elseif (trial_num==[5:8])
            end

            
    
            if size(fooof_results.peak_params,1)==1

                peak_center_table(i)=fooof_results.peak_params(1); %(1) for center (2) for power (3) for bandwidth
                peak_power_table(i)=fooof_results.peak_params(2);
                peak_bandwidth_table(i)=fooof_results.peak_params(3)/2;

                min_peak_fq(i)=peak_center_table(i)-peak_bandwidth_table(i);
                max_peak_fq(i)=peak_center_table(i)+peak_bandwidth_table(i);
                peak_power_higest(i)=peak_power_table(i);
                
            elseif size(fooof_results.peak_params,1)==2

                for peak=1:2
                    peak_center_table(peak,i)=fooof_results.peak_params(peak,1); %(1) for center (2) for power (3) for bandwidth
                    peak_power_table(peak,i)=fooof_results.peak_params(peak,2);
                    peak_bandwidth_table(peak,i)=fooof_results.peak_params(peak,3)/2;
                end

                min_peak_fq(i)=peak_center_table(1,i)-peak_bandwidth_table(1,i);
                max_peak_fq(i)=peak_center_table(2,i)+peak_bandwidth_table(2,i);
                peak_power_higest(i)=max(peak_bandwidth_table(:,i));
            else
                peak_center_table(i)=0; %(1) for center (2) for power (3) for bandwidth
                peak_power_table(i)=0;
                peak_bandwidth_table(i)=0;
                min_peak_fq(i)=peak_center_table(i)-peak_bandwidth_table(i);
                max_peak_fq(i)=peak_center_table(i)+peak_bandwidth_table(i);
                peak_power_higest(i)=peak_power_table(i);
            end

            plot_parame=[min_peak_fq(i)
                        0
                        max_peak_fq(i)-min_peak_fq(i)
                        peak_power_higest(i)];
            rectangle('Position',plot_parame','edgecolor','r','linewidth',0.8);%,'edgecolor','k','facecolor','g','linewidth',1.8

        end
        
        perodic_frequence(trial_num,:)=mean(peak_center_table(:,:),1);
        perodic_frequence_new(trial_num)=mean(perodic_frequence(trial_num,2:3),2);

        sgtitle(['use fooof get peak parameters = ' num2str(perodic_frequence_new(trial_num))])
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Fourier transfor segment  Data  to two group
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_forier
        figure
        
        for i=1:3
            dt=pause_time;% sampling interval
            N=length(plot_data(:,1));% sampling points
            t=time';% sampling time
            Fs=Fs;%Sampling frequency, which is the reciprocal of the sampling interval
            n=0:1:N-1;
            f=(Fs/N).*n;%The frequency corresponding to each point on the X axis
            x=plot_data(:,i); %signal
    
            y=fft(x);% Fourier transform to get a complex number
            Ay=abs(y);% modulo
            Ayy=Ay*2/N;% converted to actual amplitude
    
            subplot(6,3,0*3+i)
            plot(t,x);
            title('time figure')
            subplot(6,3,1*3+i)
            plot(f(1:N/2),Ayy(1:N/2))
            fshift=f;
            yshift=y;
            case_axis=3;
            %--- different axis
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
            title('frequence figure')
            sgtitle('all data')
            
            if use_fooof_get_centerAndBw==1
                f1=min_peak_fq(i);
                f2=max_peak_fq(i);
                if f1<0
                    f1=0;
                end
            else
                f1=0.5;
                f2=5;
            end
            yy=zeros(1,length(y));
            for m=0:N-1 
                if(m*(Fs/N)>f1&m*(Fs/N)&&(Fs-f2)&m*(Fs/N)<(Fs-f1));
                   yy(m+1)=0; 
                else 
                    yy(m+1)=y(m+1); 
                end
            end
    
            dd=zeros(1,length(y));
            for m=0:N-1 
                if(m*(Fs/N)<f1||m*(Fs/N)>f2);%filter out the frequency after Nyquist
                   dd(m+1)=0; 
                else 
                    dd(m+1)=y(m+1); 
                end
            end
            
            yyi=abs(yy);
            yi=ifft(yy);
    
            subplot(6,3,2*3+i)
            plot(f(1:N/2),yyi(1:N/2))
            fshift=f;
            yshift=yyi;
            %--- different axis
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
            title('frequence figure')
            subplot(6,3,3*3+i)
            plot(t,real(yi))
            title('time figure')
            sgtitle('only non-periodic part')
    
            ddi=abs(dd);
            di=ifft(dd);
            
    
            subplot(6,3,4*3+i)
            plot(f(1:N/2),ddi(1:N/2))
            fshift=f;
            yshift=ddi;
            %--- different axis
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
            title('frequence figure')
            subplot(6,3,5*3+i)
            plot(t,real(di))
            title('time figure')
            sgtitle('only periodic part')
    
            plot_data_p(:,i)=real(di);
            plot_data_np(:,i)=real(yi);
    
        end
    end

    if plot_segment==1
        figure
        subplot(2,2,[1 2])
        plot3(plot_data_p(:,1),plot_data_p(:,2),plot_data_p(:,3));box on; grid on;
        title('only non-periodic part')
        xlabel('x');ylabel('y');zlabel('z');
        subplot(223)
        plot3(plot_data(:,1),plot_data(:,2),plot_data(:,3));box on; grid on;
        title('all data')
        xlabel('x');ylabel('y');zlabel('z');
        view(68,35)
        subplot(224)
        plot3(plot_data_np(:,1),plot_data_np(:,2),plot_data_np(:,3));box on; grid on;
        title('only periodic part')
        sgtitle('see 3D data')
        xlabel('x');ylabel('y');zlabel('z');
        view(68,35)
    end

   


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save data for future training
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if save_or_not==1
        proc_data_sub.time=time;
        proc_data_sub.periodic=plot_data_p;
        proc_data_sub.non_periodic=plot_data_np;
        proc_data_sub.quant_hand=proc_data{trial_num}.Q_xyzw(1:4,:)';
        proc_data_sub.quant_target=proc_data{trial_num}.target_Qwxyz';
        proc_data_sub.angle_vel=proc_data{trial_num}.AngleVelocity(1:3,:)';
        proc_data_sub.X_dot=proc_data{trial_num}.X_dot(1:3,:)';

        proc_data_seprated{trial_num}=proc_data_sub;
    end

    
    if save_or_not==1
        clear plot_data_p plot_data_np proc_data_sub plot_data_p_fooof plot_data_np_fooof
    end

end

if use_fooof_get_centerAndBw==1
    perodic_frequence_new'
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save data for future training
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_or_not==1
    save([path_of_save exp_kind '_seprated'],'proc_data_seprated');
    !echo SAVE TO MAT DONE!


else

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot for paper
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    h=figure
    for i=2
        dt=pause_time;% sampling interval
        N=length(plot_data(:,1));% sampling points
        t=time';% sampling time
        Fs=Fs;%Sampling frequency, which is the reciprocal of the sampling interval
        n=0:1:N-1;
        f=(Fs/N).*n;%The frequency corresponding to each point on the X axis
        x=plot_data(:,i); %signal
    
        y=fft(x);% Fourier transform to get a complex number
        
        Ay=abs(y);% modulo
        Ayy=Ay*2/N;% converted to actual amplitude
    
        
        if i==1
            title('X direction')
        elseif i==2
            title('Y direction')
        elseif i==3
            title('Z direction')
        end
        
        
        
        freqs           =  f         %row vector of frequency values
        power_spectrum  =  abs(y')         %row vector of power values
        f_range         =  [0, 10]         %fitting range (Hz)
        %--- settings :    %fooof model settings, in a struct, including:
        settings.peak_width_limts= [2.0 3.0];
        settings.max_n_peaks= 1;
        settings.min_peak_height=0.2;
        settings.peak_threshold=2.0;
        settings.aperiodic_mode='fixed';
        return_model    = 1         %boolean of whether to return the FOOOF model fit, optional
        
        fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model)
    
        subplot(2,3,2)
        fooof_plot_subplot(h,fooof_results);hold on;grid on
        title('fitting results for power spectrum')
        subplot(2,3,5)
        p=plot(fooof_results.freqs,(fooof_results.fooofed_spectrum-fooof_results.ap_fit));hold on;grid on
        xlabel('Frequence');ylabel('Magnitude');
        title('bandwidth of periodic signal')
        %--- invers fuliye
        di=ifft(fooof_results.ap_fit);
        yi=ifft(fooof_results.fooofed_spectrum-fooof_results.ap_fit);
    
        plot_data_p_fooof(:,i)=real(yi);
        plot_data_np_fooof(:,i)=real(di);
    
        %--- manually set some peak
        if (trial_num==1)
            if i==1
                fooof_results.peak_params(1)=2.2;
                fooof_results.peak_params(3)=0.7;
            end
        elseif (trial_num==2)
            if i==1
                fooof_results.peak_params(1)=2.3;
                fooof_results.peak_params(3)=0.5;
            end
        elseif (trial_num==4)
            if i==1
                fooof_results.peak_params(1)=2.7;
                fooof_results.peak_params(2)=5.4;
                fooof_results.peak_params(3)=0.4;
            end
        elseif (trial_num==[5:8])
        end
    
        
    
        if size(fooof_results.peak_params,1)==1
    
            peak_center_table(i)=fooof_results.peak_params(1); %(1) for center (2) for power (3) for bandwidth
            peak_power_table(i)=fooof_results.peak_params(2);
            peak_bandwidth_table(i)=fooof_results.peak_params(3)/2;
    
            min_peak_fq(i)=peak_center_table(i)-peak_bandwidth_table(i);
            max_peak_fq(i)=peak_center_table(i)+peak_bandwidth_table(i);
            peak_power_higest(i)=peak_power_table(i);
            
        elseif size(fooof_results.peak_params,1)==2
    
            for peak=1:2
                peak_center_table(peak,i)=fooof_results.peak_params(peak,1); %(1) for center (2) for power (3) for bandwidth
                peak_power_table(peak,i)=fooof_results.peak_params(peak,2);
                peak_bandwidth_table(peak,i)=fooof_results.peak_params(peak,3)/2;
            end
    
            min_peak_fq(i)=peak_center_table(1,i)-peak_bandwidth_table(1,i);
            max_peak_fq(i)=peak_center_table(2,i)+peak_bandwidth_table(2,i);
            peak_power_higest(i)=max(peak_bandwidth_table(:,i));
        else
            peak_center_table(i)=0; %(1) for center (2) for power (3) for bandwidth
            peak_power_table(i)=0;
            peak_bandwidth_table(i)=0;
            min_peak_fq(i)=peak_center_table(i)-peak_bandwidth_table(i);
            max_peak_fq(i)=peak_center_table(i)+peak_bandwidth_table(i);
            peak_power_higest(i)=peak_power_table(i);
        end
    
        plot_parame=[min_peak_fq(i)
                    0
                    max_peak_fq(i)-min_peak_fq(i)
                    peak_power_higest(i)];
        r=rectangle('Position',plot_parame','edgecolor','r','linewidth',0.8);%,'edgecolor','k','facecolor','g','linewidth',1.8
        
        hold on;
        p_fake = plot(NaN, NaN, '-r'); 
        legend([p,p_fake],{'periodic signal','bandwidth'})
    end
    
    
    subplot(2,3,3)
    plot3(plot_data_p(:,1),plot_data_p(:,2),plot_data_p(:,3));box on; grid on;
    title('only non-periodic part')
    xlabel('x');ylabel('y');zlabel('z');
    subplot(2,3,[1 4])
    plot3(plot_data(:,1),plot_data(:,2),plot_data(:,3));box on; grid on;
    title('all data')
    xlabel('x');ylabel('y');zlabel('z');
    view(68,35)
    subplot(236)
    plot3(plot_data_np(:,1),plot_data_np(:,2),plot_data_np(:,3));box on; grid on;
    title('only periodic part')
    xlabel('x');ylabel('y');zlabel('z');
    view(68,35)
    

    Fig5_left=plot_data(:,1:3);
    Fig5_right_periodic=plot_data_p(:,1:3);
    Fig5_right_aperiodic=plot_data_np(:,1:3);
    % Saving these variables to a MAT file
    save('Fig5_left_right_Data.mat', 'Fig5_left', 'Fig5_right_periodic', 'Fig5_right_aperiodic');



    figure
    plot3(plot_data(:,1),plot_data(:,2),plot_data(:,3));box on; grid on;
    % title('nozzle motion in')
    xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');

    a1=figure
    fooof_plot_subplot(a1,fooof_results);hold on;grid on
    % title('fitting results for power spectrum')

    figure
    p=plot(fooof_results.freqs,(fooof_results.fooofed_spectrum-fooof_results.ap_fit));hold on;grid on
    xlabel('Frequence');ylabel('Magnitude');
    % title('bandwidth of periodic signal')

    figure
    plot3(plot_data_p(:,1),plot_data_p(:,2),plot_data_p(:,3));box on; grid on;
    % title('only non-periodic part')
    xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');

    figure
    plot3(plot_data_np(:,1),plot_data_np(:,2),plot_data_np(:,3));box on; grid on;
    % title('only periodic part')
    xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
    view(68,35)


end
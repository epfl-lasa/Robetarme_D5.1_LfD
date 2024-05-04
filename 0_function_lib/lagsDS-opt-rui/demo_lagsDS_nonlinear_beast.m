%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Script for GMM-based LAGS-DS Learning introduced in paper:         %
%  'Locally Active Globally Stable Dynamical Systems';                    %
% N. Figueroa and A. Billard; TRO/IJRR 2019                               %
% With this script you can load/draw 2D toy trajectories and real-world   %
% trajectories acquired via kinesthetic taching and test the different    %                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019 Learning Algorithms and Systems Laboratory,          %
% EPFL, Switzerland                                                       %
% Author:  Nadia Figueroa                                                 % 
% email:   nadia.figueroafernandez@epfl.ch                                %
% website: http://lasa.epfl.ch                                            %
%                                                                         %
% This work was supported by the EU project Cogimon H2020-ICT-23-2014.    %
%                                                                         %
% Permission is granted to copy, distribute, and/or modify this program   %
% under the terms of the GNU General Public License, version 2 or any     %
% later version published by the Free Software Foundation.                %
%                                                                         %
% This program is distributed in the hope that it will be useful, but     %
% WITHOUT ANY WARRANTY; without even the implied warranty of              %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General%
% Public License for more details                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demo Code for Non-linear LAGS-DS with multiple local regions %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 1 [OPTION 1]: Load 2D Data Drawn from GUI %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear all; clc
% %%%%%%% Choose to plot robot for simulation %%%%%%%
% with_robot   = 1;
% continuous   = 1;
% 
% if with_robot
%     % set up a simple robot and a figure that plots it
%     robot = create_simple_robot();
%     fig1 = initialize_robot_figure(robot);
%     title('Feasible Robot Workspace','Interpreter','LaTex')    
%     % Base Offset
%     base = [-1 1]';
%     % Axis limits
%     limits = [-2.5 0.5 -0.45 2.2];
% else
%     fig1 = figure('Color',[1 1 1]);
%     % Axis limits
%     limits = [-2.5 0.5 -0.45 2.2];
%     axis(limits)    
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.55, 0.2646 0.4358]);
%     grid on
% end
% 
% % Draw Reference Trajectory
% data = draw_mouse_data_on_DS(fig1, limits);
% if continuous
%     % Process Drawn Data for DS learning
%     [Data, Data_sh, att_g, x0_all, dt] = processDrawnData(data);
% else
%     % Process dis-joint data another way...    
%     att_g = [0;0];    
% end
% h_att = scatter(att_g(1),att_g(2), 150, [0 0 0],'d','Linewidth',2); hold on;
% 
% % Position/Velocity Trajectories
% M          = length(Data);
% Xi_ref     = Data(1:2,:);
% Xi_dot_ref = Data(3:end,:);
% N = 2;

% %% Compute the Mean trajectory for the Locally Active Region Estimation
% sub_sample      = 1;
% [Data_mean] = getMeanTrajectory(data, [0;0], sub_sample, dt);
% scatter(Data_mean(1,:),Data_mean(2,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 1 [OPTION 2]: Load 2D Data from Real Locomotion Datasets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;
pkg_dir = '/lagsDS-opt-rui/';
%%%%%%%%%%%%%%%%%%% Choose a Dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
choosen_dataset = 4; % 1: Demos from Gazebo Simulations (right trajectories)
                     % 2: Demos from Gazebo Simulations (left+right trajectories)
                     % 3: Demos from Real iCub 
                     % 4: Demos from iCub Narrow Passage Walking
                     % 5: Demos from iCubs Object Conveyor Belt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub_sample      = 20; % To sub-sample trajectories   
[Data, Data_sh, att_g, x0_all, dt, data, ~, ...
    ~, ~, dataset_name, box_size] = load_loco_datasets(pkg_dir, choosen_dataset, sub_sample);

% Plot Position/Velocity Trajectories %%%%%
vel_samples = 20; vel_size = 0.5; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data, att_g, vel_samples, vel_size);
axis equal;
limits = axis;
h_att = scatter(att_g(1),att_g(2), 150, [0 0 0],'d','Linewidth',2); hold on;

% Position/Velocity Trajectories
[~,M]      = size(Data);
Xi_ref     = Data(1:2,:);
Xi_dot_ref = Data(3:end,:);
N = 2;

% Compute the Mean trajectory for the Locally Active Region Estimation
% [Data_mean] = getMeanTrajectory(data, att_g, sub_sample, dt);
Data_mean=Data;

% For 2D-locomotion Datasets
switch choosen_dataset
    case 1
        %%%%% Draw Obstacle %%%%%
        rectangle('Position',[-1 1 6 1], 'FaceColor',[.85 .85 .85 0.25]); hold on;
    case 2
        %%%%% Draw Obstacle %%%%%
        rectangle('Position',[-1 1 6 1], 'FaceColor',[.85 .85 .85 0.25]); hold on;
    case 3
        %%%%% Draw Table %%%%%
        rectangle('Position',[-6.75 -2.15 0.5 0.5], 'FaceColor',[.85 .85 .85 0.25]); hold on;
        Data_mean(1:2,:) = Data_mean(1:2,:) +att_g; 
    case 4
        rectangle('Position',[-1.5 1 3 2], 'FaceColor',[.85 .85 .85 0.25]); hold on;
        rectangle('Position',[2 1 3 2], 'FaceColor',[.85 .85 .85 0.25]); hold on;     
        
    case 5
        rectangle('Position',[-0.3 -1.5 0.6 3], 'FaceColor',[.85 .85 .85 0.3]); hold on;
        rectangle('Position',[2.5 3.7 3 0.6], 'FaceColor',[.85 .85 .85 0.3]); hold on;
end
title_name = strcat('Reference Trajectories from:', dataset_name);
title(title_name,'Interpreter','LaTex','FontSize',16);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 2 (Locally Linear State-Space Paritioning): Fit GMM to Trajectory Data  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% GMM Estimation Algorithm %%%%%%%%%%%%%%%%%%%%%%
% 0: Physically-Consistent Non-Parametric (Collapsed Gibbs Sampler)
% 1: GMM-EM Model Selection via BIC
% 2: CRP-GMM (Collapsed Gibbs Sampler)
est_options = [];
est_options.type             = 0;   % GMM Estimation Alorithm Type   

% If algo 0 or 2 selected:
est_options.samplerIter      = 50;  % Maximum Sampler Iterations                                 
est_options.do_plots         = 1;   % Plot Estimation Statistics
est_options.sub_sample       = 2;   % Size of sub-sampling of trajectories
                                    % 1/2 for 2D datasets, >2/3 for real    
% Metric Hyper-parameters
est_options.estimate_l       = 1;   % '0/1' Estimate the lengthscale, if set to 1
est_options.l_sensitivity    = 2;   % lengthscale sensitivity [1-10->>100]
est_options.length_scale     = [];  

% Fit GMM to Trajectory Data
[Priors, Mu, Sigma] = fit_gmm(Xi_ref, Xi_dot_ref, est_options);
% [Priors, Mu, Sigma] = fit_gmm(Data_mean(1:2,:), Data_mean(3:4,:), est_options);
K = length(Priors);

%% Generate GMM data structure for DS learning
clear ds_gmm; ds_gmm.Mu = Mu; ds_gmm.Sigma = Sigma; ds_gmm.Priors = Priors;   

%% Align Covariance Matrices Re-estimate GMM parameters
unique_labels = unique(est_labels);
Mu_k = Mu;  Sigma_k = Sigma;
for k=1:length(unique_labels)
    cluster_points = Xi_ref(:,est_labels == unique_labels(k));
    if ~isempty(cluster_points)
        [ V_k, L_k, Mu_k(:,k) ] = my_pca( cluster_points );
        Sigma_k(:,:,k) = V_k*L_k*V_k';
    end
end
rel_dilation_fact = 0.25;
Sigma_k = adjust_Covariances(Priors, Sigma_k, 1, rel_dilation_fact);
ds_gmm.Mu = Mu_k; ds_gmm.Sigma = Sigma_k;

%%  Visualize Gaussian Components and labels on clustered trajectories 
% Extract Cluster Labels
[~, est_labels] =  my_gmm_cluster(Xi_ref, ds_gmm.Priors, ds_gmm.Mu, ds_gmm.Sigma, 'hard', []);

% Visualize Estimated Parameters
[h_gmm]  = visualizeEstimatedGMM(Xi_ref,  ds_gmm.Priors, ds_gmm.Mu, ds_gmm.Sigma, est_labels, est_options);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Step 3: Estimate Local directions/attractors and Hyper-Plane Functions  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Create Radial-Basis Function around global attractor %%%%
B_r = 10;
radius_fun = @(x)(1 - my_exp_loc_act(B_r, att_g, x));
grad_radius_fun = @(x)grad_lambda_fun(x, B_r, att_g);

%%%%    Extract local linear DS directions and attractors   %%%%
% If attractors inside local regions decrease gauss_thres
gauss_thres = 0.25; plot_local_params = 1;
[att_l, local_basis] = estimate_local_attractors_lags(Data, est_labels, ds_gmm, gauss_thres, radius_fun, att_g);
if plot_local_params
    [h_att_l, h_dirs] = plotLocalParams(att_l, local_basis, Mu);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 4: ESTIMATE CANDIDATE LYAPUNOV FUNCTION PARAMETERS  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learn Single P_g for Initialization of WSAQF
[Vxf] = learn_wsaqf_lags(Data_sh);
P_g_prior = Vxf.P;

% Learn WSAQF 
scale_prior = 2; % ==> Check that it converges.. if it doesn't, reduce scale!
att_origin = att_l-repmat(att_g,[1 K]);
att_l_mod = att_l;

[Vxf]    = learn_wsaqf_lags(Data_sh, att_origin, P_g_prior*scale_prior);
P_g      = Vxf.P(:,:,1);
P_l      = Vxf.P(:,:,2:end);

%%% If any local attractor is equivalent to the global re-parametrize %%%
att_diff = att_l-repmat(att_g,[1 K]);
equal_g = find(any(att_diff)==0);
if equal_g ~= 0
    for ii=1:length(equal_g)
        P_g_ = P_g + P_l(:,:,equal_g(ii));
        P_l(:,:,equal_g(ii)) = P_g_;
    end
end
%% %%% Plot learned Lyapunov Function %%%%%
if N == 2
    contour = 0; % 0: surf, 1: contour
    clear lyap_fun grad_lyap 
    
    % Lyapunov function
    lyap_fun = @(x)lyapunov_function_combined(x, att_g, att_l_mod, 1, P_g, P_l, ds_gmm);
    % Gradient of Lyapunov function
    grad_lyap = @(x)gradient_lyapunov_combined(x, att_g, att_l_mod, P_g, P_l);
            
    title_string = {'$V(\xi) = (\xi-\xi_g^*)^TP_g(\xi-\xi_g^*) + \sum_{k=1}^K\beta^k((\xi-\xi_g^*)^TP_l^k(\xi-\xi_k^*))^2$'};
%     if exist('h_lyap','var');  delete(h_lyap);     end
%     if exist('hd_lyap','var');     delete(hd_lyap);     end
    [h_lyap] = plot_lyap_fct(lyap_fun, contour, limits,  title_string, 0);
    [hd_lyap] = scatter(Data(1,:),Data(2,:),10,[1 1 0],'filled'); hold on;            
    [h_grad_lyap] = plot_gradient_fct(grad_lyap, limits,  '$V(\xi)$ and $\nabla_{\xi}V(\xi)$ Function');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Step 5: Learn Global DS Parameters     %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choose Global DS type %%%
globalDS_type = 1; % 0: linear-DS
                   % 1: nonlinear-DS, LPV-DS learned with WSAQF Lyapunov Function 
%%% Choose Lyapunov Candidate type %%%                   
lyapType      = 1; % 0: P-QLF as in CoRL paper
                   % 1: WSAQF as in LAGS paper

clear global_ds
switch globalDS_type
    case 0 % Parametrize an agnostic global DS
    %%%%%%%%  Linear DS f_g(xi) = (Axi + b) %%%%%%%%
    
    globalDS_name = 'Linear Global DS';
    case 1 
    %%%%%%%%  LPV-DS f_g(xi) = sum_{k=1}^{K}\gamma_k(xi)(A_kxi + b_k) %%%%%%%%                
        switch lyapType
            case 0
                %%%%%%% Learn LPV-DS with P-QLF Constraint %%%%%
                % Type of constraints/optimization
                constr_type = 2; init_cvx    = 0;                                
                [A_g, b_g, ~] = optimize_lpv_ds_from_data(Data, att_g, constr_type, ds_gmm, P_g, init_cvx);                
                global_ds = @(x) lpv_ds(x, ds_gmm, A_g, b_g);                
                globalDS_name = 'Global LPV-DS (P-QLF)';
                
            case 1
                %%%%%%% Learn LPV-DS with WSAQF Constraint %%%%%
                eps_scale = 0.1; % ==> "If Constraints don't converge increase/decrease"
                enforce_g = 1;  % ==> "Enforce A+A' If local constraints not met
                
                % Slow down the velocitues for global DS
                Data_         = Data; scale = 5;
                Data_(3:4,:)  = scale*Data_(3:4,:);
                [A_g, b_g, ~] = optimize_globalDS_lags(Data_, att_g, 3, ds_gmm, P_g, P_l, att_l_mod, eps_scale, enforce_g, equal_g);                                
                global_ds = @(x) (1/scale)*lpv_ds(x, ds_gmm, A_g, b_g); 
                globalDS_name = 'Global LPV-DS (WSAQF)';                                                              
        end        
end

%% %%%  Plot Resulting DS  %%%%%
% Fill in plotting options
ds_plot_options = [];
ds_plot_options.sim_traj  = 1;            % To simulate trajectories from x0_all
ds_plot_options.x0_all    = x0_all;       % Intial Points
% ds_plot_options.x0_all    = x0_all_new;       % Intial Points
ds_plot_options.init_type = 'ellipsoid';  % For 3D DS, to initialize streamlines
                                          % 'ellipsoid' or 'cube'  
ds_plot_options.nb_points = 30;           % No of streamlines to plot (3D)
ds_plot_options.plot_vol  = 1;            % Plot volume of initial points (3D)
% ds_plot_options.limits    = limits;

[hd, hs, hr, x_sim] = visualizeEstimatedDS(Xi_ref, global_ds, ds_plot_options);
limits = axis;
title(globalDS_name, 'Interpreter','LaTex','FontSize',20)
% [hatt_rob] = scatter(Data_mean(1,:),Data_mean(2,:), 150, [0 1 0],'s','Linewidth',2); hold on;

%% %%%%%%%%%%%%   Export DS parameters to Mat/Txt/Yaml files  %%%%%%%%%%%%%%%%%%%
DS_name = 'iCub-Object-Conveyor';
save_lpvDS_to_Mat(DS_name, pkg_dir, ds_gmm, A_g, b_g, att_g, x0_all, dt, P_g, 2, est_options)

%% Save LPV-DS parameters to text files
DS_name = 'iCub-Object-Conveyor';
save_lpvDS_to_txt(DS_name, pkg_dir,  ds_gmm, (1/scale)*A_g, att_g)

%% Save LPV-DS parameters to yaml file
DS_name = 'iCub-Object-Conveyor';
% To use the rest of the code you need a matlab yaml convertor
% you can get it from here: http://vision.is.tohoku.ac.jp/~kyamagu/software/yaml/
save_lpvDS_to_Yaml(DS_name, pkg_dir,  ds_gmm, (1/scale)*A_g, att_g, x0_all, dt)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      [OPTIONAL]: Evaluate Global DS Accuracy/Stability   %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% Compare Velocities from Demonstration vs DS %%%%%
% Compute RMSE on training data
rmse = mean(rmse_error(global_ds, Xi_ref, Xi_dot_ref));
fprintf('Global-DS has prediction RMSE on training set: %d \n', rmse);

% Compute e_dot on training data
edot = mean(edot_error(global_ds, Xi_ref, Xi_dot_ref));
fprintf('Global-DS has e_dot on training set: %d \n', edot);

if exist('h_vel_comp','var'); delete(h_vel_comp); end
[h_vel_comp] = visualizeEstimatedVelocities(Data, global_ds);
title('Real vs. Estimated Velocities w/Global-DS', 'Interpreter', 'LaTex', 'FontSize', 15)

%% %%%  Plot Lyapunov function Derivative %%%%%
if N == 2
    clear lyap_der
    contour = 0; %0 :surface, 1:contour 
    switch lyapType
        case 0            
            % Lyapunov function
            lyap_fun = @(x)lyapunov_function_PQLF(x, att_g, P_g);
            % Plots
            title_string = {'$V(\xi) = (\xi-\xi^*)^TP(\xi-\xi^*)$'};
            h_lyap     = plot_lyap_fct(lyap_fun, contour, limits,  title_string, 0);
            [hddd] = scatter(Data(1,:),Data(2,:),10,[1 1 0],'filled'); hold on;
    
            % Derivative of Lyapunov function (gradV*f(x))
            lyap_der = @(x)lyapunov_derivative_PQLF(x, att_g, P_g, global_ds);
            title_string_der = {'Lyapunov Function Derivative $\dot{V}(\xi)$'};
        case 1
            % Derivative of Lyapunov function (gradV*f(x))
            lyap_der = @(x)lyapunov_combined_derivative(x, att_g, att_l_mod, global_ds, 1, P_g, P_l);                                  
            title_string_der = {'Lyapunov Function Derivative $\dot{V}_{DD-WSAQF}(\xi)$'};
    end
%     if exist('h_lyap_der','var'); delete(h_lyap_der); end
%     if exist('hd_lyap_der','var'); delete(hd_lyap_der); end
    [h_lyap_der] = plot_lyap_fct(lyap_der, contour, limits,  title_string_der, 1);
    [hd_lyap_der] = scatter(Data(1,:),Data(2,:),10,[1 1 0],'filled'); hold on;            
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      [OPTIONAL]: Simulate Global DS on Robot with Passive-DS Ctrl   %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up a simple robot and a figure that plots it
robot = create_simple_robot();
fig1 = initialize_robot_figure(robot);
title('Feasible Robot Workspace','Interpreter','LaTex')
% Base Offset
base = [-1 1]';
% Axis limits
limits = [-2.5 0.75 -0.45 2.2];
axis(limits)
if exist('hds_rob','var');  delete(hds_rob);  end
if exist('hdata_rob','var');delete(hdata_rob);end
if exist('hatt_rob','var');delete(hatt_rob);end
[hds_rob] = plot_ds_model(fig1, global_ds, [0;0], limits,'medium'); hold on;
[hdata_rob] = scatter(Data(1,:),Data(2,:),10,[1 0 0],'filled'); hold on;            
[hatt_rob] = scatter(att_g(1),att_g(2), 150, [0 0 0],'d','Linewidth',2); hold on;
title(globalDS_name, 'Interpreter','LaTex','FontSize',20)
    xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
    ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);
    
%% Simulate Passive DS Controller function
dt = 0.01;
show_DK = 0;
if show_DK
    struct_stiff = [];
    struct_stiff.DS_type  = 'global'; % Options: 'global', 'lags'
    struct_stiff.gmm      = ds_gmm;
    struct_stiff.A_g      = A_g;
    struct_stiff.basis    = 'D'; % Options: D (using the same basis as D) or I (using the world basis)
    struct_stiff.L        = [0.5 0; 0 0.5]; % Eigenvalues for Daming matrix
    struct_stiff.mod_step = 50; % Visualize D/K matrices every 'mod_step'     
    
    % Select this one if we want to see the Damping and Stiffness matrices
    simulate_passiveDS(fig1, robot, base, global_ds, att_g, dt, struct_stiff);
    
else
    simulate_passiveDS(fig1, robot, base, global_ds, att_g, dt);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 6: Create Activation Function (i.e. select locallly active regions)     %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Global DS with Locally Linear Partitions
clc;
ds_plot_options.sim_traj  = 0;            
[hd, hs, hr, x_sim, fig_ds] = visualizeEstimatedDS(Data_mean(1:2,:), global_ds, ds_plot_options);
title('Choose regions to locally activate:', 'Interpreter','LaTex','FontSize',20);
[ h_gmm, h_ctr, h_txt ] = plotLLRegions(ds_gmm.Mu, ds_gmm.Sigma);

% Compute Hyper-plane functions per local model
w = zeros(2,K); breadth_mod = 20;
h_functor = cell(1,K); grad_h_functor = cell(1,K);
lambda_functor  = cell(1,K); Data_k = [];
for k=1:K
    % Create Hyper-Plane functions
    % Insert gradient of lyapunov function in one of these functions!
    w(:,k)                 = -local_basis(:,1,k);
    h_functor{k}           = @(x)hyper_plane(x,w(:,k),att_l(:,k));
    grad_h_functor{k}      = @(x)grad_hyper_plane(x,w(:,k),h_functor{k});
    lambda_functor{k}      = @(x)(1-my_exp_loc_act(breadth_mod, att_l(:,k), x));
%     lambda_functor{k}      = @(x)lambda_mod_fun(x, breadth_mod, att_l(:,k), grad_h_functor{k}, grad_lyap);
    grad_lambda_functor{k} = @(x)grad_lambda_fun(x, breadth_mod, att_l(:,k));    
end

%% Choose Locally Linear Regions that you wanto to activate

% Choose the regions that you want to locally activate
choosen_active = 1:K;
% choosen_active = 2:K;
% choosen_active = [4 5];
% choosen_active = [2 3 4];
% choosen_active = [4 8 5 9 7];
% choosen_active = [2 5 3];
% choosen_active = [1 2 3 4 8 5];

% Construct Dataset for those regions
[~, est_labels_mean] =  my_gmm_cluster(Data_mean(1:2,:), ds_gmm.Priors, ds_gmm.Mu, ds_gmm.Sigma, 'hard', []);

%% Construct Activation Function
act_type = 1; % 0: GMM-based
              % 1: GPR-based  

% Radial Basis Function centered at the global attractor
B_r = 5;
radius_fun = @(x)(1 - my_exp_loc_act(B_r, att_g, x));
grad_radius_fun = @(x)grad_lambda_fun(x, B_r, att_g);

% Construct GMM or GPR activation functions
switch act_type
    case 0 % With GMM                
        
        % Hyper-params
        peak_scale = 0.98;
        
        % Parametrize Priors, Mu, Sigma for GMM-based Activation Function
        K_act       = length(choosen_active);
        Priors_act  = (1/K_act)*ones(1,K_act);                
        Mu_act      = ds_gmm.Mu(:,choosen_active);
        Sigma_act   = ds_gmm.Sigma(:,:,choosen_active)*0.9;
        
        % Generate GMM data structure for Mapping Function
        clear gmr_fun activ_fun alpha_fun
        gmm_fun    = @(x)ml_gmm_pdf(x, Priors_act, Mu_act, Sigma_act);
        min_peak   = min(gmm_fun(Mu_act))*peak_scale;
        activ_fun  = @(x) 1 - min(min_peak,gmm_fun(x))./min_peak;
        alpha_fun  = @(x)( (1-radius_fun(x))'.*activ_fun(x)' + radius_fun(x)');         
                
        % Compute gradient       
        grad_gmm_fun   = @(x)grad_gmm_pdf(x, ds_gmm);
        grad_alpha_fun = @(x)gradient_alpha_fun(x,radius_fun, grad_radius_fun, gmm_fun, grad_gmm_fun, 'gmm');
        activation_name = 'GMM';
        
    case 1 % With GPR      
        
        % Hyper-params
        rbf_var = 0.35;
        
        % Construct Data for GPR
        Data_act = []; 
        for k_=1:length(choosen_active)
            Data_act       = [Data_act Data_mean(1:2,est_labels_mean==choosen_active(k_))];
        end
        if length(Data_act)>500
            sub_sample_gpr = floor(length(Data_act)/200);
        else
            sub_sample_gpr= 1;
        end            
        Data_act = Data_act(:,1:sub_sample_gpr:end);
        X_act = Data_act; y_act = ones(1,length(X_act));
        
        % Parametrize GPR
        epsilon = 0.0001; % assuming output variance = 1
        model.X_train   = X_act';   model.y_train   = y_act';
        clear gpr_fun activ_fun alpha_funSwipe
        gpr_fun    = @(x) my_gpr(x',[],model,epsilon,rbf_var);
        activ_fun  = @(x) 1 - gpr_fun(x);                    
        alpha_fun  = @(x)max(0,(1-radius_fun(x))'.*activ_fun(x)' + radius_fun(x)');        
        
        % Compute gradient       
        clear grad_gpr_fun grad_alpha_fun
        grad_gpr_fun   = @(x)gradient_gpr(x, model, epsilon, rbf_var);
        grad_alpha_fun = @(x)gradient_alpha_fun(x,radius_fun, grad_radius_fun, gpr_fun, grad_gpr_fun, 'gpr');        
        activation_name = 'GPR';        
end

%% Plot Activation function
figure('Color',[1 1 1])
[h_act] = plot_mixing_fct_2d(limits, alpha_fun); hold on;
[hds_rob] = plot_ds_model(h_act, global_ds, [0;0], limits,'medium'); hold on;
[hdata_act] = scatter(Data(1,:),Data(2,:),10,[0 0 0],'filled'); hold on;            
[hatt_act] = scatter(att_g(1),att_g(2), 150, [0 0 0],'d','Linewidth',2); hold on;
[ h_gmm, h_ctr, h_txt ] = plotLLRegions(ds_gmm.Mu, ds_gmm.Sigma, choosen_active);
[h_att_l, h_dirs] = plotLocalParams(att_l, local_basis, Mu, choosen_active);
xlabel('$\xi_1$','Interpreter', 'LaTex','FontSize',15)
ylabel('$\xi_2$','Interpreter', 'LaTex','FontSize',15)
title_name = strcat(activation_name, '-based Activation function');
title(title_name, 'Interpreter','LaTex','FontSize',20);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Step 7: Learn Local DS Parameters      %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Option 1: Estimate Local Dynamics from Data by searching for optimal Kappa via parallel convex optimization problems %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_l_k = zeros(N,N,K); b_l_k = zeros(N,K);
A_d_k = zeros(N,N,K); b_d_k = zeros(N,K);
show_plots = 0; kappas = zeros(1,K);
clc;
for k=1:K     
    if any(find(choosen_active==k))
        fprintf('\n****************************************************');
        fprintf('\n******* Starting %d-th local-DS Optimization *******',k);
        fprintf('\n****************************************************');
        Data_k = Data_mean(:,est_labels_mean==k);
        sub_sample = 10;
        Data_k = Data_k(:,1:sub_sample:end);
        
        % Initial evaluation points
        kappas_sample = [1, 15, 30];
        DS_g = (k == equal_g);
        i = 1;
        if DS_g
            kappas_sample = [1, 10, 20]            
        end
        
        [A_l_i(:,:,i), b_l_i(:,i), A_d_i(:,:,i), b_d_i(:,i), num_grid_violations] = estimate_localDS_search(Data_k, kappas_sample(i), 1/scale*A_g(:,:,k),  att_g, att_l(:,k), P_g, P_l(:,:,k), ...
                    local_basis(:,:,k), alpha_fun, h_functor{k}, grad_lyap,grad_h_functor{k}, activ_fun, ds_gmm.Mu(:,k), ds_gmm.Sigma(:,:,k), lambda_functor{k}, DS_g);
        viols_sample = zeros(1,3);
        
        
        viols_sample(1,i) = num_grid_violations;
        if viols_sample(1,i) > 0
            % If DS is not stable at lowest kappa, estimate local-DS with full constraints
            %%%%%%%%%%%%%%%%%%%%%%  LOCAL DS OPTIMIZATION OPTIONS %%%%%%%%%%%%%%%%%%%%%
            % Construct variables and function handles for stability constraints
            clear stability_vars
            stability_vars.solver          = 'fmincon'; % options: 'baron' or 'fmincon'
            stability_vars.grad_h_fun      = grad_h_functor{k};
            stability_vars.add_constr      = 1; % 0: no stability constraints
            
            % 1: Adding sampled stability constraints
            % Type of constraint to evaluate
            stability_vars.constraint_type = 'full';   % options: 'full/matrix/hessian'
            stability_vars.epsilon         = -150;     % 'small' number for f_Q < -eps
            stability_vars.do_plots        = 0;        % plot current lyapunov constr.
            stability_vars.init_samples    = 100;      % Initial num/boundary samples
            stability_vars.iter_samples    = 250;      % Initial num/boundary samples
            
            % Function handles for contraint evaluation
            stability_vars.alpha_fun       = alpha_fun;
            stability_vars.grad_alpha_fun  = grad_alpha_fun;
            stability_vars.activ_fun       = activ_fun;
            stability_vars.h_fun           = h_functor{k};
            stability_vars.lambda_fun      = lambda_functor{k};
            stability_vars.P_l             = P_l;
            stability_vars.P_g             = P_g;
            
            % Variable for different type constraint types
            if strcmp(stability_vars.constraint_type,'full')
                stability_vars.grad_lyap_fun = grad_lyap;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf('\n=========== Estimating %d-th Local Dynamics ===========\n',k);
            %%%%%%%%%% To compute outer contour for Chi-samples %%%%%%%%%%
            Norm = ml_gaussPDF(Mu(:,k),Mu(:,k),Sigma(:,:,k))*0.75;
            [A_l, b_l, A_d, b_d] = estimate_localDS_multi(Data_k, 1/scale*A_g(:,:,k), att_g, 1, att_l(:,k), local_basis(:,:,k), ds_gmm.Mu(:,k), ds_gmm.Sigma(:,:,k), Norm, limits, stability_vars);
            Lambda_l = eig(A_l);
            kappa = max(abs(Lambda_l))/min(abs(Lambda_l));
        else
            for i=2:3
                [A_l_i(:,:,i), b_l_i(:,i), A_d_i(:,:,i), b_d_i(:,i), num_grid_violations] = estimate_localDS_search(Data_k, kappas_sample(i), 1/scale*A_g(:,:,k),  att_g, att_l(:,k), P_g, P_l(:,:,k), ...
                    local_basis(:,:,k), alpha_fun, h_functor{k}, grad_lyap,grad_h_functor{k}, activ_fun, ds_gmm.Mu(:,k), ds_gmm.Sigma(:,:,k),lambda_functor{k}, DS_g);
                viols_sample(1,i) = num_grid_violations;
            end
                        
            % If DS is stable at max kappa
            if viols_sample(3) == 0
                A_l = A_l_i(:,:,3); A_d = A_d_i(:,:,3);
                b_l = b_l_i(:,3);   b_d = b_d_i(:,3);
                kappa = kappas_sample(1,3);
                
            elseif viols_sample(2) == 0
                %%%%%%%%%%%%%%%%%%%%%%  LOCAL DS OPTIMIZATION OPTIONS %%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%% WITH PARAMETER SEARCH STRATEGY %%%%%%%%%%%%%%%%%%%%%%
                
                fprintf('\n***** Searching for 25 < kappa < 50! *******');
                new_kappa = (kappas_sample(3)-kappas_sample(2))/2 + kappas_sample(2);
                % Initialize search
                [A_l, b_l, A_d, b_d, num_grid_violations] = estimate_localDS_search(Data_k, new_kappa, 1/scale*A_g(:,:,k),  att_g, att_l(:,k), P_g, P_l(:,:,k), ...
                    local_basis(:,:,k), alpha_fun, h_functor{k}, grad_lyap,grad_h_functor{k}, activ_fun, ds_gmm.Mu(:,k), ds_gmm.Sigma(:,:,k),lambda_functor{k}, DS_g);
                clear est_params
                est_params.A_l = A_l; est_params.b_l = b_l;
                est_params.A_d = A_d; est_params.b_d = b_d;
                est_params.kappa_min          = kappas_sample(2);
                est_params.kappa_max          = new_kappa;
                est_params.current_violations = num_grid_violations;
                [A_l, b_l, A_d, b_d, num_grid_violations,new_kappa] = search_localDS_Range(est_params, Data_k, 1/scale*A_g(:,:,k), att_g, att_l(:,k), P_g, P_l(:,:,k), ...
                    local_basis(:,:,k), alpha_fun, h_functor{k}, grad_lyap,grad_h_functor{k}, activ_fun, ds_gmm.Mu(:,k), ds_gmm.Sigma(:,:,k),lambda_functor{k});
                kappa = new_kappa;
                
            else
                %%%%%%%%%%%%%%%%%%%%%%  LOCAL DS OPTIMIZATION OPTIONS %%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%% WITH PARAMETER SEARCH STRATEGY %%%%%%%%%%%%%%%%%%%%%%
                fprintf('\n***** Searching for 1 < kappa < 25! *******');
                new_kappa = kappas_sample(2)/2;
                % Initialize search
                [A_l, b_l, A_d, b_d, num_grid_violations] = estimate_localDS_search(Data_k, new_kappa, 1/scale*A_g(:,:,k),  att_g, att_l(:,k), P_g, P_l(:,:,k), ...
                    local_basis(:,:,k), alpha_fun, h_functor{k}, grad_lyap,grad_h_functor{k}, activ_fun, ds_gmm.Mu(:,k), ds_gmm.Sigma(:,:,k),lambda_functor{k}, DS_g);
                clear est_params
                est_params.A_l = A_l; est_params.b_l = b_l;
                est_params.A_d = A_d; est_params.b_d = b_d;
                est_params.kappa_min          = 1;
                est_params.kappa_max          = new_kappa;
                est_params.current_violations = num_grid_violations;
                [A_l, b_l, A_d, b_d, num_grid_violations,new_kappa] = search_localDS_Range(est_params, Data_k, 1/scale*A_g(:,:,k), att_g, att_l(:,k), P_g, P_l(:,:,k), ...
                    local_basis(:,:,k), alpha_fun, h_functor{k}, grad_lyap,grad_h_functor{k}, activ_fun, ds_gmm.Mu(:,k), ds_gmm.Sigma(:,:,k),lambda_functor{k});
                kappa = new_kappa;
            end           
        end
        fprintf('\n=========== DONE %d-th DS with kappa=%2.4f ===========\n',k, kappa);
        if show_plots
            linearDS_k = @(x) linearDS(x, A_l, b_l);
            fig_k = figure('Color',[1 1 1]);
            [hs] = plot_ds_model(fig_k, linearDS_k, [0 0]', limits,'medium'); hold on;
            [hdata_rob] = scatter(Data_k(1,:),Data_k(2,:),10,[0 0 0],'filled'); hold on;
            box on
            grid on
            xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
            ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);
        end
        
        % Fill in full matrices
        A_l_k(:,:,k) = A_l; b_l_k(:,k) = b_l;
        A_d_k(:,:,k) = A_d; b_d_k(:,k) = b_d;
        kappas(1,k) = kappa;
    end    
end
fprintf('\n******* Finished local-DS Optimization *******\n');
kappas_search = kappas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Option 2:  Constrained Local-Dynamics Estimation, estimating \kappa_k with full constraint via NLP %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappas = zeros(1,K);
A_l_k = zeros(N,N,K); b_l_k = zeros(N,K);
A_d_k = zeros(N,N,K); b_d_k = zeros(N,K);
show_plots = 0;
clc;
for k=1:K    
    if any(find(choosen_active==k))
        Data_k = Data_mean(:,est_labels_mean==k);
        
        %%%%%%%%%%%%%%%%%%%%%%  LOCAL DS OPTIMIZATION OPTIONS %%%%%%%%%%%%%%%%%%%%%
        % Construct variables and function handles for stability constraints
        clear stability_vars
        stability_vars.solver          = 'baron'; % options: 'baron' or 'fmincon'
        stability_vars.grad_h_fun      = grad_h_functor{k};
        stability_vars.add_constr      = 1; % 0: no stability constraints
        
        % 1: Adding sampled stability constraints
        % Type of constraint to evaluate
        stability_vars.constraint_type = 'full';   % options: 'full/matrix/hessian'
        stability_vars.epsilon         = -150;     % small number for f_Q < -eps
        stability_vars.do_plots        = 0;        % plot current lyapunov constr.
        stability_vars.init_samples    = 100;      % Initial num/boundary samples
        stability_vars.iter_samples    = 250;      % Initial num/boundary samples
        
        % Function handles for contraint evaluation
        stability_vars.alpha_fun       = alpha_fun;
        stability_vars.grad_alpha_fun  = grad_alpha_fun;
        stability_vars.activ_fun       = activ_fun;
        stability_vars.h_fun           = h_functor{k};
        stability_vars.lambda_fun      = lambda_functor{k};
        stability_vars.P_l             = P_l;
        stability_vars.P_g             = P_g;
        
        % Variable for different type constraint types
        if strcmp(stability_vars.constraint_type,'full')
            stability_vars.grad_lyap_fun = grad_lyap;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fprintf('\n=========== Estimating %d-th Local Dynamics ===========\n',k);        
        %%%%%%%%%% To compute outer contour for Chi-samples %%%%%%%%%%
        Norm = ml_gaussPDF(Mu(:,k),Mu(:,k),Sigma(:,:,k))*0.75;        
        [A_l, b_l, A_d, b_d] = estimate_localDS_multi(Data_k, A_g(:,:,k), att_g, 1, att_l(:,k), local_basis(:,:,k), ds_gmm.Mu(:,k), ds_gmm.Sigma(:,:,k), Norm, limits, stability_vars);                                
        Lambda_l = eig(A_l);
        kappas(1,k) = max(abs(Lambda_l))/min(abs(Lambda_l));                  
                       
        fprintf('\n=========== DONE %d-th DS with kappa=%2.4f ===========\n',k, kappas(1,k));
        
        % Fill in full matrices
        A_l_k(:,:,k) = A_l; b_l_k(:,k) = b_l;
        A_d_k(:,:,k) = A_d; b_d_k(:,k) = b_d; 
        
    end
end
fprintf('\n******* Finished local-DS Optimization *******\n');
kappas_nlp = kappas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Option 3: Manually set the convergence rates   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds_types = ones(1,K);
tracking_factor = 15*ones(1,K); % Relative ratio between eigenvalues
A_l_k = zeros(N,N,K); b_l_k = zeros(N,K);
A_d_k = zeros(N,N,K); b_d_k = zeros(N,K);
for k=1:K
    if any(find(choosen_active==k))
    Data_k = Data_mean(:,est_labels_mean==k);
    sub_sample = 1;
    Data_k = Data_k(:,1:sub_sample:end);
    
    % Compute lamda's from reference trajectories
%     [A_l_k(:,:,k), ~, ~, ~] = estimate_localDS_known_gamma(Data_k, 1/scale*A_g(:,:,k),  att_g, att_l(:,k), ds_types(k), tracking_factor(k),w(:,k),P_g, P_l(:,:,k), local_basis(:,:,k));
%     b_l_k(:,k) = -A_l_k(:,:,k)*att_l(:,k);
        
    % Manual A_l using the rayleigh quotient
    kappa = tracking_factor(k);
    A_g_k = 1/scale*A_g(:,:,k);    
    R_gk = (local_basis(:,1,k)'*A_g_k*local_basis(:,1,k))/(local_basis(:,1,k)'*local_basis(:,1,k));    
   
    A_l_k(:,:,k) = local_basis(:,:,k)*[0.75*R_gk 0;0 0.75*kappa*R_gk]*local_basis(:,:,k)'
    b_l_k(:,k) = -A_l_k(:,:,k)*att_l(:,k);
    
    % Check Eigenvalues
    Lambda_l = eig(A_l_k(:,:,k));
    
    % Construct the Deflective terms
    A_d_k(:,:,k) = -max(Lambda_l)*eye(N);
    b_d_k(:,k) = -A_d_k(:,:,k)*att_l(:,k);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Create Function for DS and PLot                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modulation = 1;
lags_ds = @(x) lags_ds_nonlinear(x, alpha_fun, ds_gmm, 1/scale*A_g, 1/scale*b_g, A_l_k, b_l_k, A_d_k, b_d_k, h_functor, lambda_functor, grad_h_functor, att_g, att_l, modulation, scale);

%%%%%  Plot Resulting DS  %%%%%
fig_lagds = figure('Color',[1 1 1]);
[h_act] = plot_mixing_fct_2d(limits, alpha_fun); hold on;
[hds] = plot_ds_model(fig_lagds, lags_ds, [0;0], limits,'medium'); hold on;
[hdata_ds] = scatter(Data_mean(1,:),Data_mean(2,:),10,[0 0 1],'filled'); hold on;            
% [hdata_ds] = scatter(Data(1,:),Data(2,:),10,[1 0 0],'filled'); hold on;

%%%%%  Generate Simulations from DS  %%%%%
% Example 4
% x0_all_new = [x0_all [-1.31;0.25] [-0.74;0.235] [-0.6396;1.533] [-1.218;1.283] [-0.8142;-0.4516] [-0.6414;1.689]];
% x0_all_new = [x0_all [-1.31;0.25] [-0.74;0.235] [-0.6396;1.533] [-1.218;1.283] [-0.8142;-0.4516] [-0.6414;1.689]];

% Example 5
% x0_all_new = [x0_all [-1.31;0.25]   [-0.77 ;0]];
if true
x0_all_new = [x0_all];
opt_sim = [];
opt_sim.dt    = 0.05;
opt_sim.i_max = 100000;
opt_sim.tol   = 0.005;
opt_sim.plot  = 0;
[x_sim, ~]    = Simulation(x0_all_new ,[],@(x)(2*lags_ds(x)), opt_sim);
[hdata_rob] = scatter(x_sim(1,:),x_sim(2,:),20,[0 0 0],'filled'); hold on;            
[hatt_rob] = scatter(x0_all_new(1,:),x0_all_new(2,:), 150, [0 1 0],'s','Linewidth',2); hold on;
end
[hatt_rob] = scatter(att_g(1),att_g(2), 150, [0 0 0],'d','Linewidth',2); hold on;
title('Non-Linear LAGS-DS with Multi-Active Region', 'Interpreter','LaTex','FontSize',20)
xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);


%% %%%  Plot Lyapunov function %%%%%
if N == 2
    clear lyap_der
    contour = 0; %0 :surface, 1:contour lagsDS
    switch lyapType
        case 0            
            % Lyapunov function
            lyap_fun = @(x)lyapunov_function_PQLF(x, att_g, P_g);
            % Plots
            title_string = {'$V(\xi) = (\xi-\xi^*)^TP(\xi-\xi^*)$'};
            h_lyap     = plot_lyap_fct(lyap_fun, contour, limits,  title_string, 0);
            [hddd] = scatter(Data(1,:),Data(2,:),10,[1 1 0],'filled'); hold on;
    
            % Derivative of Lyapunov function (gradV*f(x))
            lyap_der = @(x)lyapunov_derivative_PQLF(x, att_g, P_g, lags_ds);
            title_string_der = {'Lyapunov Function Derivative $\dot{V}(\xi)$'};
        case 1
            % Derivative of Lyapunov function (gradV*f(x))
            lyap_der = @(x)lyapunov_combined_derivative(x, att_g, att_l_mod, lags_ds, 1, P_g, P_l);                       
            
            title_string_der = {'Lyapunov Function Derivative $\dot{V}_{DD-WSAQF}(\xi)$'};
    end
    if exist('h_lyap_der','var'); delete(h_lyap_der); end
    if exist('hd_lyap_der','var'); delete(hd_lyap_der); end
    [h_lyap_der] = plot_lyap_fct(lyap_der, contour, limits,  title_string_der, 1);
    [hd_lyap_der] = scatter(Data(1,:),Data(2,:),10,[1 1 0],'filled'); hold on;            
end

%% Save LAGS-DS parameters to text files
% DS_name = 'iCub-Narrow-Passage-LAGS';
DS_name = 'iCub-Object-Conveyor-LAGS';
save_lagsDS_to_txt(DS_name, pkg_dir,  ds_gmm, (1/scale)*A_g, att_g, A_l_k, A_d_k, att_l, w, breadth_mod*ones(1,K), scale, B_r)

% Write Model to txt file for C++ Class
clear gpr_model
gpr_model.D            = N;
gpr_model.N_train      = length(Data_act);
gpr_model.x_train      = X_act;
gpr_model.y_train      = y_act;
gpr_model.length_scale = rbf_var;
gpr_model.sigma_n      = epsilon;
gpr_model.sigma_f      = 1;

model_dir = strcat(pkg_dir,'/models/',DS_name, '/');
gpr_filename = strcat(model_dir,'GPR_model.txt');
writeGPRModel(gpr_model, gpr_filename)


%% Save LAGS-DS parameters to yaml file
% To use the rest of the code you need a matlab yaml convertor
% you can get it from here: http://vision.is.tohoku.ac.jp/~kyamagu/software/yaml/
save_lagsDS_to_Yaml(DS_name, pkg_dir,  ds_gmm, (1/scale)*A_g, att_g, x0_all, dt, A_l_k, A_d_k, att_l, w, breadth_mod*ones(1,K), scale, B_r, gpr_filename)

%% Optional save reference trajectories with computed velocities for C++ class testing

% Create function handle for local component estimation only
lags_ds_local  = @(x) lags_ds_nonlinear_local(x, alpha_fun, ds_gmm, A_g, b_g, A_l_k, b_l_k, A_d_k, b_d_k, h_functor, lambda_functor, grad_h_functor, att_g, att_l, 1, scale);

% Randomly sample points within the range of the DS
N_samples = 200;
Xmin_values   = [-0.015;-0.015];
Xrange_values = [4;4];
Data_grid = Xmin_values(:,ones(1,N_samples)) + rand(N_samples,N)'.*(Xrange_values(:, ones(1,N_samples)));
% scatter(Data_grid(1,:), Data_grid(2,:));
Data_grid(N+1:N+2,:) = zeros(N, length(Data_grid));

% Simulate velocities from same reference trajectory
Data_test = [Data Data_grid];
% Data_test = [Data ];
xd_dot_global = global_ds(Data_test(1:N,:));    
xd_dot_local  = lags_ds_local(Data_test(1:N,:));
xd_alpha      = alpha_fun(Data_test(1:N,:));
xd_dot_lags   = lags_ds(Data_test(1:N,:));

% Writing Data
model_dir = strcat(pkg_dir,'/models/',DS_name, '/');

% Writing xi_ref
dlmwrite(strcat(model_dir,'Data'), Data_test, 'newline','unix','Delimiter',' ','precision','%.8f');
% Writing xi_dot_g
dlmwrite(strcat(model_dir,'xi_dot_g'), xd_dot_global, 'newline','unix','Delimiter',' ','precision','%.8f');
% Writing xi_dot_l
dlmwrite(strcat(model_dir,'xi_dot_l'), xd_dot_local, 'newline','unix','Delimiter',' ','precision','%.8f');
% Writing xi_alpha
dlmwrite(strcat(model_dir,'xi_alpha'), xd_alpha, 'newline','unix','Delimiter',' ','precision','%.8f');
% Writing xi_dot
dlmwrite(strcat(model_dir,'xi_dot_lags'), xd_dot_lags, 'newline','unix','Delimiter',' ','precision','%.8f');


%% Write Test Data to txt file for C++ Class [only needed for GPR]
x_test = Data_test(1:N,:);
y_test = gpr_fun(x_test);
model_dir = strcat(pkg_dir,'/models/',DS_name, '/');
gpr_testfilename = strcat(model_dir,'GPR_data.txt');
writeGPRTestData(x_test, y_test, gpr_testfilename)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      [OPTIONAL]: Simulate LAGS DS on Robot with Passive-DS Ctrl   %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up a simple robot and a figure that plots it
robot = create_simple_robot();
fig1 = initialize_robot_figure(robot);
title('Feasible Robot Workspace','Interpreter','LaTex')
% Base Offset
base = [-1 1]';
% Axis limits
% limits = [-2.5 0.75 -0.45 2.2];
limits = [-2.5 0.75 -0.45 1.5];
axis(limits)
if exist('hds_rob','var');  delete(hds_rob);  end
if exist('hdata_rob','var');delete(hdata_rob);end
if exist('hatt_rob','var');delete(hatt_rob);end
% if exist('h_act','var');delete(h_act);end
[h_act] = plot_mixing_fct_2d(limits, alpha_fun); hold on;
[hds_rob] = plot_ds_model(fig1, lags_ds, [0;0], limits,'medium'); hold on;
[hdata_rob] = scatter(Data(1,:),Data(2,:),10,[1 0 0],'filled'); hold on;            
[hdata_rob] = scatter(Data_mean(1,:),Data_mean(2,:),10,[0 0 1],'filled'); hold on;            
[hatt_rob] = scatter(att_g(1),att_g(2), 150, [0 0 0],'d','Linewidth',2); hold on;
title('Non-Linear LAGS-DS', 'Interpreter','LaTex','FontSize',20)
    xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
    ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);
% axis tight    
%% Simulate Passive DS Controller function
dt = 0.01;
show_DK = 1;
if show_DK
    struct_stiff = [];
    struct_stiff.DS_type             = 'lags'; % Options: 'global', 'lags'
    struct_stiff.gmm                 = ds_gmm;
    struct_stiff.A_g                 = A_g;   
    struct_stiff.A_l                 = A_l_k;
    struct_stiff.A_d                 = A_d_k;
    struct_stiff.att_k               = att_l;    
    struct_stiff.alpha_fun           = alpha_fun;
    struct_stiff.grad_alpha_fun      = grad_alpha_fun;
    struct_stiff.h_functor           = h_functor;
    struct_stiff.grad_h_functor      = grad_h_functor;
    struct_stiff.lambda_functor      = lambda_functor;
    struct_stiff.grad_lambda_functor = grad_lambda_functor;    
    struct_stiff.basis               = 'D';  % Options: D (using the same basis as D) or I (using the world basis)
    struct_stiff.L                   = [1 0; 0 2]; % Eigenvalues for Damping matrix
    struct_stiff.mod_step            = 10; % Visualize D/K matrices every 'mod_step'     
    
    % Select this one if we want to see the Damping and Stiffness matrices
    simulate_passiveDS(fig1, robot, base, lags_ds, att_g, dt, struct_stiff);
    
else
    simulate_passiveDS(fig1, robot, base, lags_ds, att_g, dt);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      [OPTIONAL]: Evaluate Global DS Accuracy/Stability   %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% Compare Velocities from Demonstration vs DS %%%%%
% Compute RMSE on training data
rmse = mean(rmse_error(lags_ds, Xi_ref, Xi_dot_ref));
fprintf('Global-DS has prediction RMSE on training set: %d \n', rmse);

% Compute e_dot on training data
edot = mean(edot_error(lags_ds, Xi_ref, Xi_dot_ref));
fprintf('Global-DS has e_dot on training set: %d \n', edot);

if exist('h_vel_comp','var'); delete(h_vel_comp); end
[h_vel_comp] = visualizeEstimatedVelocities(Data, lags_ds);
title('Real vs. Estimated Velocities w/LAGS-DS', 'Interpreter', 'LaTex', 'FontSize', 15)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    [OPTIONAL] 3D Visualization of any function      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize a function
eval_fun       = @(x)alpha_fun(x);
eval_fun       = @(x)radius_fun(x);
eval_fun       = @(x)lambda_functor{1}(x)
% eval_fun = @(x)vecnorm(x-att_g);
% modulation = 1;
% norm_ds_loc  = @(x) lags_ds_nonlinear_local(x, alpha_fun, ds_gmm, A_g, b_g, A_l_k, b_l_k, A_d_k, b_d_k, h_functor, lambda_functor, grad_h_functor, att_g, att_l, modulation, scale);
% norm_ds_glob = @(x)vecnorm(x-att_g);
% gamma = 1;
% eval_fun = @(x)(norm_ds_loc(x)-gamma*norm_ds_glob(x));

% eval_fun = lambda_functor{1};
h_gamma      = plot_lyap_fct(eval_fun, 0, limits,  {'$\alpha(\xi)$ Function'}, 1);
[hdata_act] = scatter(Data(1,:),Data(2,:),10,[1 0 0],'filled'); hold on;            
% Visualize its gradient
% grad_eval_fun  = @(x)grad_alpha_fun(x);
% h_grad_gamma = plot_gradient_fct(grad_eval_fun, limits,  '$r$ and $\nabla_{\xi}r$ Function');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Script for GMM-based LAGS-DS Learning introduced in paper:         %
%  'Locally Active Globally Stable Dynamical Systems';                    %
% N. Figueroa and A. Billard; IJRR 2019                               %
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 1: Load 2D Data from Real Locomotion Datasets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

pkg_dir_load     = './0_human_demo';
pkg_dir_save     = './0_function_lib/lagsDS-opt-rui';

DS_name         ='2023-05-23-snake-shape';
chosen_dataset  = 1; 
sub_sample      = 1; % '>2' for real 3D Datasets, '1' for 2D toy datasets
nb_trajectories = 3; % For real 3D data only
[Data, Data_sh, att, x0_all, data, dt] = load_dataset_DS(pkg_dir_load, chosen_dataset, sub_sample, nb_trajectories);

att_g = att;

%%---- velocity some time is not good
velocity_not_good=1;
dt=0.01;
if velocity_not_good==1
    for i=1:3
        Data(i+3,:)=gradient(Data(i,:),dt);
    end
end
if velocity_not_good==1
    for i=1:3
        Data_sh(i+3,:)=gradient(Data_sh(i,:),dt);
    end
end

%%%%% Plot Position/Velocity Trajectories %%%%%
vel_samples = 15; vel_size = 0.5; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data, att, vel_samples, vel_size);
axis equal;
limits = axis;
h_att = scatter3(att(1),att(2),att(3), 150, [0 0 0],'d','Linewidth',2); hold on;

% Extract Position and Velocities
M          = size(Data,1)/2;    
Xi_ref     = Data(1:M,:);
Xi_dot_ref = Data(M+1:end,:);   


% Compute the Mean trajectory for the Locally Active Region Estimation
sub_sample = 2;
[Data_mean] = getMeanTrajectory(data, att_g, sub_sample, dt);

scatter3(Data_mean(1,:),Data_mean(2,:),Data_mean(3,:), 'k*');

limits = axis;

title_name = strcat('Reference Trajectories');
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
est_options.sub_sample       = 1;   % Size of sub-sampling of trajectories
                                    % 1/2 for 2D datasets, >2/3 for real    
% Metric Hyper-parameters
est_options.estimate_l       = 1;   % '0/1' Estimate the lengthscale, if set to 1
est_options.l_sensitivity    = 10;   % lengthscale sensitivity [1-10->>100]
est_options.length_scale     = [];  

% Fit GMM to Trajectory Data
% [Priors, Mu, Sigma] = fit_gmm(Xi_ref, Xi_dot_ref, est_options);
[Priors, Mu, Sigma] = fit_gmm(Xi_ref, Xi_dot_ref, est_options);
K = length(Priors);

%% Generate GMM data structure for DS learning
clear ds_gmm; 
ds_gmm.Mu = Mu; ds_gmm.Sigma = Sigma; ds_gmm.Priors = Priors;   

%% Align Covariance Matrices Re-estimate GMM parameters
% Extract Cluster Labels
[~, est_labels] =  my_gmm_cluster(Xi_ref, ds_gmm.Priors, ds_gmm.Mu, ds_gmm.Sigma, 'hard', []);
unique_labels = unique(est_labels);
Mu_k = Mu;  Sigma_k = Sigma;
for k=1:length(unique_labels)
    cluster_points = Xi_ref(:,est_labels == unique_labels(k));
    if ~isempty(cluster_points)
        [ V_k, L_k, Mu_k(:,k) ] = my_pca( cluster_points );
        Sigma_k(:,:,k) = V_k*L_k*V_k';
    end
end
rel_dilation_fact = 0.35;
Sigma_k = adjust_Covariances(Priors, Sigma_k, 1, rel_dilation_fact);
ds_gmm.Mu = Mu_k; ds_gmm.Sigma = Sigma_k;

%%  Visualize Gaussian Components and labels on clustered trajectories 
% Visualize Estimated Parameters
[h_gmm]  = visualizeEstimatedGMM(Xi_ref,  ds_gmm.Priors, ds_gmm.Mu, ds_gmm.Sigma, est_labels, est_options);
hold on;
scatter3(Data_mean(1,:),Data_mean(2,:),Data_mean(3,:), 'k*')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Step 3: Estimate Local directions/attractors and Hyper-Plane Functions  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Create Radial-Basis Function around global attractor %%%%
% Real iCub Datasets
B_r = 50;
% KUKA with LASA DATASET
% B_r = 200;
radius_fun = @(x)(1 - my_exp_loc_act(B_r, att_g, x));
grad_radius_fun = @(x)grad_lambda_fun(x, B_r, att_g);

clc;
%%%%    Extract local linear DS directions and attractors   %%%%
% If attractors inside local regions decrease gauss_thres
gauss_thres = 0.1; 
% gauss_thres = 20; 

plot_local_params = 0;
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
scale_prior = 1; % ==> Check that it converges.. if it doesn't, reduce scale!
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

    contour = 0; % 0: surf, 1: contour
    clear lyap_fun grad_lyap 
    
    % Lyapunov function
    lyap_fun = @(x)lyapunov_function_combined(x, att_g, att_l_mod, 1, P_g, P_l, ds_gmm);
    
    % Gradient of Lyapunov function
    grad_lyap = @(x)gradient_lyapunov_combined(x, att_g, att_l_mod, P_g, P_l);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Step 5: Learn Global DS Parameters     %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choose Global DS type %%%
globalDS_type = 1; % 0: linear-DS
                   % 1: nonlinear-DS, LPV-DS learned with WSAQF Lyapunov Function 
%%% Choose Lyapunov Candidate type %%%                   
lyapType      = 0; % 0: P-QLF as in CoRL paper
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
                %%%%%%---------
                scale = 10;
                %%%%%%%%%%% always cause some error
                
                [A_g, b_g, ~] = optimize_lpv_ds_from_data(Data, att_g, constr_type, ds_gmm, P_g_prior, init_cvx);                
                global_ds = @(x) lpv_ds(x, ds_gmm, A_g, b_g);                
                globalDS_name = 'Global LPV-DS (P-QLF)';
                
            case 1
                %%%%%%% Learn LPV-DS with WSAQF Constraint %%%%%
                eps_scale = 1; % ==> "If Constraints don't converge increase/decrease"
                enforce_g = 1;  % ==> "Enforce A+A' If local constraints not met
                
                % Slow down the velocities for global DS
                Data_         = Data; scale = 10;
                Data_(4:6,:)  = scale*Data_(4:6,:);
                [A_g, b_g, ~] = optimize_globalDS_lags(Data_, att_g, 3, ds_gmm, P_g, P_l, att_l_mod, eps_scale, enforce_g, equal_g);                                
                global_ds = @(x) (1/scale)*lpv_ds(x, ds_gmm, A_g, b_g); 
                globalDS_name = 'Global LPV-DS (WSAQF)';                                                              
        end        
end


%% %%%  Plot Resulting DS  %%%%%
% Fill in plotting options
ds_plot_options = [];
ds_plot_options.sim_traj  = 1;            % To simulate trajectories from x0_all
if ds_plot_options.sim_traj == 1
    ds_plot_options.x0_all    = x0_all;       % Intial Points
    ds_plot_options.init_type = 'ellipsoid';  % For 3D DS, to initialize streamlines
    % 'ellipsoid' or 'cube'
    ds_plot_options.nb_points = 30;           % No of streamlines to plot (3D)
    ds_plot_options.plot_vol  = 1;            % Plot volume of initial points (3D)
else
    ds_plot_options.x0_all    = [];
end
ds_plot_options.limits    = limits;

[hd, hs, hr, x_sim] = visualizeEstimatedDS(Xi_ref, global_ds, ds_plot_options);
limits = axis;
title(globalDS_name, 'Interpreter','LaTex','FontSize',20)
[hatt_rob] = scatter3(Data_mean(1,:),Data_mean(2,:),Data_mean(3,:), 50, [1 0 0],'s','Linewidth',2); hold on;
h_att = scatter3(att_g(1),att_g(2),att_g(3),150, [0 0 0],'d','Linewidth',2); hold on;
text(att_g(1),att_g(2),att_g(3),'$\mathbf{\xi}^*_g$','Interpreter', 'LaTex','FontSize',15); hold on;

%% plot for 2D : Plot YZ coordinates of DS
yz = [2 3];
ds_plot_options = [];
ds_plot_options.sim_traj   = 1;      % To simulate trajectories from x0_all
ds_plot_options.x0_all     = x0_all; % Initial Points for reproductions
ds_plot_options.dimensions = yz;
ds_plot_options.attractor  = att;
[hd, hs, hr, x_sim] = visualizeEstimatedDS(Xi_ref, global_ds, ds_plot_options);
title('YZ-slice of learned DS', 'Interpreter','LaTex','FontSize',20)


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
if M == 2
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
    [h_lyap_der] = plot_lyap_fct(lyap_der, contour, limits,  title_string_der, 1);
    [hd_lyap_der] = scatter(Data(1,:),Data(2,:),10,[1 1 0],'filled'); hold on;            
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 6: Create Activation Function (i.e. select locallly active regions)     %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Global DS with Locally Linear Partitions
clc;
ds_plot_options.sim_traj  = 0;            
% [hd, hs, hr, x_sim, fig_ds] = visualizeEstimatedDS(Xi_ref, global_ds, ds_plot_options);
[hd, hs, hr, x_sim] = visualizeEstimatedDS(Xi_ref, global_ds, ds_plot_options);
title('Choose regions to locally activate:', 'Interpreter','LaTex','FontSize',20);

%%%%%%%---- break inside this function to find out which GMM we need
%%%%%%%active---%%%
[ h_gmm, h_ctr, h_txt ] = plotLLRegions(ds_gmm.Mu, ds_gmm.Sigma);

% Compute Hyper-plane functions per local model
w = zeros(M,K); breadth_mod = 100;
h_functor = cell(1,K); grad_h_functor = cell(1,K);
lambda_functor  = cell(1,K); Data_k = [];
for k=1:K
% for k=1
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

%--- Choose the regions that you want to locally activate
choosen_active = 1:K
choosen_active = [1 2 6 7 5 9];

%--- Construct Dataset for those regions
[~, est_labels_mean] =  my_gmm_cluster(Xi_ref, ds_gmm.Priors, ds_gmm.Mu, ds_gmm.Sigma, 'hard', []);

%% Construct Activation Function
act_type = 1; % 0: GMM-based
              % 1: GPR-based  

% Radial Basis Function centered at the global attractor
% Real iCub Datasets
B_r = 50;
% real robot [lasa-dataset]
B_r_g = 1; 
radius_fun = @(x)(1 - my_exp_loc_act(B_r_g, att_g, x));
grad_radius_fun = @(x)grad_lambda_fun(x, B_r_g, att_g);

% Construct Data used for activation function
Data_act = [];  Data_act_vel = [];
for k_=1:length(choosen_active)
    Data_act       = [Data_act Xi_ref(:,est_labels_mean==choosen_active(k_))];
    Data_act_vel    = [Data_act_vel Xi_dot_ref(:,est_labels_mean==choosen_active(k_))];
end
Data_act = Data_act(:,10:end-5);
Data_act_vel = Data_act_vel(:,10:end-5);


% Construct GMM or GPR activation functions
switch act_type
    case 0 % With GMM                
        
        % Hyper-params
        peak_scale = 0.85;
        
        % Parametrize Priors, Mu, Sigma for GMM-based Activation Function
        K_act       = length(choosen_active);
        Priors_act  = (1/K_act)*ones(1,K_act);                
        Mu_act      = ds_gmm.Mu(:,choosen_active);
        Sigma_act   = ds_gmm.Sigma(:,:,choosen_active)*1.5;
        
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
        rbf_var = 0.35; % WR: adjust the size of local activity regen
       
        max_gpr_samples = 100;
        if length(Data_act)>max_gpr_samples
            sub_sample_gpr = floor(length(Data_act)/max_gpr_samples);
        else
            sub_sample_gpr= 1;
        end            
        Data_act = Data_act(:,1:sub_sample_gpr:end);
        Data_act_vel = Data_act_vel(:,1:sub_sample_gpr:end);
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

% Plot Activation function
if M==2
figure('Color',[1 1 1])
[h_act] = plot_mixing_fct_3d(limits, alpha_fun); hold on;
[hds_rob] = plot_ds_model(h_act, global_ds, [0;0], limits,'medium'); hold on;
[hdata_act] = scatter(Data(1,:),Data(2,:),10,[1 0 0],'filled'); hold on;            
[hatt_act] = scatter(att_g(1),att_g(2), 150, [0 0 0],'d','Linewidth',2); hold on;

[hdata_act] = scatter(Data_act(1,:),Data_act(2,:),50,[0 0 0],'filled'); hold on;            
[ h_gmm, h_ctr, h_txt ] = plotLLRegions(ds_gmm.Mu, ds_gmm.Sigma, choosen_active);
[h_att_l, h_dirs] = plotLocalParams(att_l, local_basis, Mu, choosen_active);

xlabel('$\xi_1$','Interpreter', 'LaTex','FontSize',15)
ylabel('$\xi_2$','Interpreter', 'LaTex','FontSize',15)
title_name = strcat(activation_name, '-based Activation function');
title(title_name, 'Interpreter','LaTex','FontSize',20);
end

%% Visualize in 3D
if M==2
eval_fun       = @(x)alpha_fun(x);
h_gamma      = plot_lyap_fct(eval_fun, 0, limits,  {'$\alpha(\xi)$ Function'}, 1);
[hdata_act] = scatter(Data(1,:),Data(2,:),10,[1 0 0],'filled'); hold on;            

zlabel('$\alpha(\xi)$', 'Interpreter','LaTex','FontSize',20);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save learned result for parameter change
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save workspace incase lost data
% save([pkg_dir_save '/models/' DS_name]);
load([pkg_dir_save '/models/' DS_name]);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Step 7: Learn Local DS Parameters      %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Option 3: Manually set the convergence rates   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds_types = ones(1,K);


tracking_factor = ones(1,K)*100
% tracking_factor=[500 300 100 1500 300 15000 15000 1500 300 6000]*0.3;
tracking_factor(choosen_active) = 300
tracking_factor([6 7]) = 1500


A_l_k = zeros(M,M,K); b_l_k = zeros(M,K);
A_d_k = zeros(M,M,K); b_d_k = zeros(M,K);
for k=1:K
    if any(find(choosen_active==k))
    Data_k = Data(:,est_labels_mean==k);
    sub_sample = 1;
    Data_k = Data_k(:,1:sub_sample:end);
    
    % Manual A_l using the rayleigh quotient
    kappa = tracking_factor(k);
    A_g_k = 1/scale*A_g(:,:,k);    
    R_gk = (local_basis(:,1,k)'*0.5*(A_g_k+A_g_k')*local_basis(:,1,k))/(local_basis(:,1,k)'*local_basis(:,1,k));    
    if R_gk > 0 
        R_gk = min(eig(0.5*(A_g_k+A_g_k')));
    end
   
    A_l_k(:,:,k) = local_basis(:,:,k)*[0.75*R_gk 0 0;0 0.75*kappa*R_gk 0; 0 0 0.75*kappa*R_gk]*local_basis(:,:,k)';
    b_l_k(:,k) = -A_l_k(:,:,k)*att_l(:,k);
    
    % Check Eigenvalues
    Lambda_l = eig(A_l_k(:,:,k));
    
    % Construct the Deflective terms
    A_d_k(:,:,k) = -max(Lambda_l)*eye(M);
    b_d_k(:,k) = -A_d_k(:,:,k)*att_l(:,k);
    end
end
tracking_factor

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Step 8: Visualize LAGS-DS      %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Create Function for DS and PLot                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modulation = 0; do_scale = 1;
% lags_ds = @(x) lags_ds_nonlinear(x, alpha_fun, d_gmm, 1/scale*A_g, 1/scale*b_g, A_l_k, b_l_k, A_d_k, b_d_k, h_functor, lambda_functor, grad_h_functor, att_g, att_l, modulation, scale, do_scale);
lags_ds = @(x) lags_ds_nonlinear(x, alpha_fun, ds_gmm, 1/scale*A_g, 1/scale*b_g, A_l_k, b_l_k, A_d_k, b_d_k, h_functor, lambda_functor, grad_h_functor, att_g, att_l, modulation, scale);

%%%%%  Plot Resulting DS  %%%%%
fig_lagds = figure('Color',[1 1 1]);
if M==2
[h_act] = plot_mixing_fct_3d(limits, alpha_fun); hold on;
end
tic;
[hds] = plot_ds_model_3D(fig_lagds, lags_ds, att_g, limits,x0_all,'low'); hold on;
toc;

% [hdata_ds] = scatter(Data_mean(1,:),Data_mean(2,:),50,[1 0 0],'filled'); hold on;            
[hdata_ds] = scatter3(Data(1,:),Data(2,:),Data(3,:),10,[1 0 0],'filled'); hold on;

%%%%%  Generate Simulations from DS  %%%%%
if false
x0_all_new = [x0_all Data(1:3,1)];
opt_sim = [];
opt_sim.dt    = 0.15;
opt_sim.i_max = 10000;
opt_sim.tol   = 0.001;
opt_sim.plot  = 0;
[x_sim, ~]    = Simulation(x0_all_new ,[],@(x)(lags_ds(x)), opt_sim);
[hdata_rob] = scatter3(x_sim(1,:),x_sim(2,:),x_sim(3,:),20,[0 0 0],'filled'); hold on;            
[hatt_rob] = scatter3(x0_all_new(1,:),x0_all_new(2,:),x0_all_new(3,:), 150, [0 1 0],'s','Linewidth',2); hold on;
end
[hatt_rob] = scatter3(att_g(1),att_g(2),att_g(3), 150, [0 0 0],'d','Linewidth',2); hold on;
title('Non-Linear LAGS-DS with Multi-Active Region', 'Interpreter','LaTex','FontSize',20)
xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);

h_att = scatter3(att_g(1),att_g(2),att_g(3), 150, [0 0 0],'d','Linewidth',2); hold on;
% legend({'$\xi^{ref}$'},'Interpreter', 'LaTex','FontSize',15)
text(att_g(1),att_g(2),att_g(3),'$\mathbf{\xi}^*_g$','Interpreter', 'LaTex','FontSize',15); hold on;

axis tight
box on
%% %%%  Plot Lyapunov function %%%%%
if M == 2
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      [OPTIONAL]: WRITE LAGS PARAMETERS TO YAML/TXT FORMAT FOR ROBOT   %%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save LAGS-DS parameters to text files
DS_name = 'wr-try-lags1';
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

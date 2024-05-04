%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rui Wu 2022.07.24
%   read human demos data and segment into peridico and nonpredico part
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all; close all;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set path
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_of_data='2ndSessionRobetarmeRecording';

if isunix
    %-----  path for linux
    path_of_load = ['./0_human_demo/' folder_of_data '/processed/'];
    path_of_plot=['./0_figure/' folder_of_data '/'];
    path_of_save = ['./0_human_demo/' folder_of_data '/trained/'];
else
    path_of_load = ['.\0_human_demo\' folder_of_data '\processed\'];
    path_of_plot=['.\0_figure\' folder_of_data '\'];
    path_of_save = ['.\0_human_demo\' folder_of_data '\trained\'];
end

status = mkdir(path_of_plot);   %path for save figure
status = mkdir(path_of_save); 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read data
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data_kind='subject1_all_perdico_data';
% data_kind='subject2_all_perdico_data';
data_kind='subject3_all_perdico_data';

load([path_of_load data_kind])

X = traj.data;
time = mean(diff(traj.time));
type = 4;
smoothing = 10;
T = [];
plot3(X(:,1),X(:,2),X(:,3));
%% % Now we prepare the data and convert it to spherical coordinates:
% Get data with preferred structure
[Xdata,Xvel,Rdata,Rvel,dt,T,N,m,begin] = prepareData(type,X,time,smoothing,T);
%% Plot data
%%
plotData(Xdata,Xvel,Rdata,Rvel,T,m,'all');
legend('Data trajectory','Location','northeast');

Pose_downsample=[Xvel(:,2:3) Xdata(:,2:3)];
vector_param.vel_samples=0.1;
vector_param.vel_size=3;
vector_param.color=[0 1 0];
plot_vector_color_2D(Pose_downsample',vector_param)

plot_vel([Xdata,Xvel]',1,1);
%% Find $\omega$ limit set
% Initialize and perform Expectation Maximization on a Gaussian Mixture Model 
% with 1 or 2 Gaussians to find the model of the limit set:
%%
% j = input('Choose number of Gaussian models (defaults to 1):\n');
% if j~=1 && j~=2
%     j = 1
% end
j = 1

[Priors, Mu, Sigma] = EM_init_kmeans([Xdata';Xvel'], j);
[Priors, Mu, Sigma] = EM([Xdata';Xvel'], Priors, Mu, Sigma);
%% 
% Plot the cluster found as the $\omega$ limit set:

figure; grid on; hold on; view(3);
title('Original data with \omega limit set found');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
if j == 2
    plotGMM(Mu(1:3,1), Sigma(1:3,1:3,1), [1 0.5 0.5], 4);
    plotGMM(Mu(1:3,2), Sigma(1:3,1:3,2), [1 0.7 0.5], 4);
    plot3(Xdata(:,1), Xdata(:,2), Xdata(:,3), 'r');
    legend('1st Gaussian model','2nd Gaussian model','Data trajectories','Location','northeast');
else
    plotGMM(Mu(1:3,:), Sigma(1:3,1:3,:), [1 0.5 0.5], 4);
    plot3(Xdata(:,1), Xdata(:,2), Xdata(:,3), 'r');
    legend('Gaussian model','Data trajectories','Location','northeast');
end
%% If wanted, perfom dimensionality reduction
% Decide if you want to perform dimensionality reduction:
%%
% prompt = 'Enter 1 to perform dimensionality reduction and 0 to skip it:\n';
% dimred = input(prompt);
dimred=1
%% 
% If wanted, perform dimensionality reduction (PCA) and plot projected data:

if (j == 2)     % if there are more clusters, choose which one to perform eigenvalue decomposition on
    prompt = 'Choose which gaussian to keep (1 for 1st, 2 for 2nd):\n';
    k = input(prompt);
else
    k = 1;
end
if (dimred == 1)
    % Get rotation matrix from PCA performed on covariance matrix of gaussian:
    [Rrot,~] = eig(Sigma(1:N,1:N,k));
    Rrot = Rrot(:,N:-1:1);
    
    % Plot projected (rotated) data -- takes a lot of time! If you trust it, set "if false"
    if true
        figure; hold on; grid on;
        subplot(1,2,1); hold on; grid on;
        title('Original data'); xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
        subplot(1,2,2); hold on; grid on;
        title('Projected data'); xlabel('e_1'); ylabel('e_2'); zlabel('e_3');
        
        for i = 1:sum(T)
            Xplot = (Rrot \ (Xdata(i,:)' - Mu(1:N,k)))';
            if N == 2
                subplot(1,2,1); plot(Xdata(i,1), Xdata(i,2), 'r.'); grid on; hold on;
                subplot(1,2,2); plot(Xplot(1), Xplot(2), 'r.'); grid on; hold on;
            else
                subplot(1,2,1); view(3); plot3(Xdata(i,1), Xdata(i,2), Xdata(i,3), 'r.'); grid on; hold on;axis equal;
                subplot(1,2,2); view(3); plot3(Xplot(1), Xplot(2), Xplot(3), 'r.'); grid on; hold on;axis equal;
            end
        end
    end
else
    % If you are not performing dimensionality reduction
    Rrot = eye(N);
end

% Find Euler angle (or rotation angle if N = 2)
if(N == 3)
    theta0 = rotm2eul(Rrot)
elseif(N == 2)
    theta0 = acos(Rrot(1,1))
end

% Get rotated data and save original data
Xdata_ = Xdata;
%% Optimization
% Select initial values of paramaters and optimize:
%%
Xdata = Xdata_;
Xdata = (Rrot \ (Xdata' - Mu(1:N,k)))';

initial_parameters = [];
% Set initial rho0 to variance of Gaussian model
initial_parameters.rho0 = 3*mean(diag(Sigma(1:N,1:N,1)));
% Set initial values of other parameters
% initial_parameters.M = 2;
% Select which parameters to optimize in first objective (M, rho0, (always 0 for R,) a, x0):
initial_parameters.first = [1 1 0 1 0];
% Select which parameters to optimize in second objective (M, (always 0 for rho0,) R, a, x0):
initial_parameters.second = [1 0 1 1 0];        % e.g. [1 1 1 1 0] ecludes x0

%%%%%% OPTIMIZATION FUNCTION: %%%%%%
[params] = optimizePars(initial_parameters,Xdata,dt,begin,10);
%% Get parameters learned
%%
rho0 = params.rho0;
M = params.M;
R = params.R;
a = params.a;

%---- Wr: define by hand
params.rho0 = 0.05;
params.M = 100.4;
% params.R = -5;
% params.a = [1 1.6125 1];


% Add the mean of the Gaussian model to x0
if isfield(initial_parameters,'x0') && isempty(initial_parameters.x0)
    x0 = -Mu(1:N,k)';
else
    x0 = (Rrot * params.x0' - Mu(1:N,k))';
end
params.x0 = x0;
params.Rrot = Rrot;
disp(params);
%% Plot learned dynamics with original data and new trajectories
%%
% Get unrotated data
Xdata = Xdata_;

% Plot original data
f=figure; 
% title('Dynamical flow in cartesian coordinates'); 
hold on; grid on;
if N == 2
    plot(Xdata(:,1),Xdata(:,2),'r.'); hold on;
else
    view(3);
    plot3(Xdata(1:T(1),1),Xdata(1:T(1),2),Xdata(1:T(1),3),'r--'); hold on;
    for i = 2:m
        plot3(Xdata(sum(T(1:i-1))+1:sum(T(1:i)),1),...
            Xdata(sum(T(1:i-1))+1:sum(T(1:i)),2),...
            Xdata(sum(T(1:i-1))+1:sum(T(1:i)),3),'r--'); hold on;
    end
end

% Plot streamlines / arrows to show dynamics
ylabel('x_2'); yl = ylim;
xlabel('x_1'); xl = xlim;
if N > 2
    zlabel('x_3'); zl = zlim;
    ngrid = 7;
end

if N == 2
    [Xs,Ys] = meshgrid(linspace(xl(1),xl(2),20),linspace(yl(1),yl(2),20));
    X_plot = [Xs(:), Ys(:)];
else
    [Xs,Ys,Zs] = meshgrid(linspace(xl(1),xl(2),ngrid),...
        linspace(yl(1),yl(2),ngrid),linspace(zl(1),zl(2),ngrid));
    X_plot = [Xs(:), Ys(:), Zs(:)];
end
Y = zeros(size(X_plot));
for i= 1:size(X_plot,1)
    [r,dr] = DS(X_plot(i,:),params);
    Y(i,1:N) = sph2cartvelocities(r,dr);
end

if N == 2
    streamslice(Xs,Ys,reshape(Y(:,1),20,20),...
        reshape(Y(:,2),20,20),'method','cubic');
else
    quiver3(X_plot(:,1),X_plot(:,2),X_plot(:,3),Y(:,1),Y(:,2),Y(:,3),'color','blue');
end

% Test dynamics for T time steps
% X0 = Xdata(1,:);
X0 = Xdata(1,:)+[ 0.001 0.01 0.02];
% X0 = [0.01 0.01 0.01];
X_s = []; Xvel_s = [];
for j = 1:size(X0,1)
    X = X0(j,:);
    for i = 1:T(1)
        X_prev = X;
        %%%%%% FUNCTION TO CALCULATE POLAR/SPHERICAL VELOCITIES: %%%%%%
        [r,dr] = DS(X,params);
        %%%%%% INTEGRATE THE DS TO GET NEXT POLAR/SPHERICAL POSITION: %%%%%%
        next_r = r + dr*dt;
        % Get next cartesian position
        X = (Rrot*(hyper2cart(next_r)./a)')' - x0;
        X_s = [X_s; X];
        Xvel_s = [Xvel_s; sph2cartvelocities(r,dr)];
        if N == 2
            plot(X(1),X(2),'k.'); hold on; grid on;
        else
            plot3(X(1),X(2),X(3),'k.'); hold on; grid on;
        end
    end
end
axis equal;box on;
view(74,16)

plot_options.x_name='x [mm]';
plot_options.y_name='y [m]';
plot_options.z_name='z [m]';
% plot_options.lengend_name={'viture attractor','real trajectory'}
% plot_options.title_name='learned LAGSDS';
plot_options.savefilename='circleDSmodela';

set_figure_for_two_clum(f,plot_options)

%% Evaluation
%% %%--- not work for low vision matlab
normalized = 1;
rmserr = RMSErr(Xdata,Xvel,params,dt,normalized);
disp('RMSE (normalized):');
disp(rmserr);

coserr = cosSim(Xdata,Xvel,params,dt);
disp('Cosine similarity:');
disp(coserr);


function [Data, Data_sh, att, x0_all, data, dt] = load_dataset_DS(pkg_dir, dataset, sub_sample, nb_trajectories)

dataset_name = [];
switch dataset
    case 1
        dataset_name = '2023-05-23-snake-shape';
end

if isempty(sub_sample)
   sub_sample = 2; 
end

% For the messy-snake dataset which is already at the origin
data_ = load(strcat(pkg_dir,'/2ndSessionRobetarmeRecording/processed/',dataset_name));
data_ = data_.Hand_dataset;
N = length(data_);    
data = []; 
traj = randsample(N, nb_trajectories)';
for l=1:nb_trajectories
    data{l} = data_{traj(l)}(:,1:sub_sample:end);        
end
[Data, Data_sh, att, x0_all, dt, data] = processDataStructure(data);

end
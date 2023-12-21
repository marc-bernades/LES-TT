%% Load filt ensemble-average datasets
try
    Data_loaded = load(["Data/" + name_file_out_filt_dataset + ".mat"]);
catch
    Data_loaded = load(["Data/Filtered_DNS_ensemble_average.mat"]);
end
Data        = Data_loaded.Data_output;

% Unpack all valriables
for ii = 1:length(Data(1,:))
    value = Data{2:end,ii};
%     value = reshape(value,[num_points_x,num_points_y,num_points_z]);
    assignin('base',Data{1,ii}, value)
end

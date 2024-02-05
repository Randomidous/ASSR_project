%% create warpmats

evLatency_left = evLatency_ms_left./2;      % from miliseconds to samples
evLatency_right = evLatency_ms_right./2;    % from miliseconds to samples
newLatency = round(maxLatency_ms./2);       % from miliseconds to samples
newLatency_ms = maxLatency_ms;

warpmats_left = cell(length(EEG_39_left.epoch),1);
for i = 1:length(EEG_39_left.epoch)
    warpmats_left{i} = timewarp(evLatency_left(i,:), newLatency);
end
warpmats_right = cell(length(EEG_39_right.epoch),1);
for i = 1:length(EEG_39_right.epoch)
    warpmats_right{i} = timewarp(evLatency_right(i,:), newLatency);
end

%% warp data

% warp left turns
warped_data_39_left = cell(length(EEG_39_left.epoch),1);     % initialize cell array
warped_data_41_left = cell(length(EEG_41_left.epoch),1);
for i = 1:length(EEG_39_left.epoch)
    data_to_warp_39_left = EEG_39_left.data(:,:,i);
    data_to_warp_39_left(:,evLatency_left(i,end)+1:end) = [];
    warped_data_39_left{i} = warpmats_left{i} * data_to_warp_39_left';
    data_to_warp_41_left = EEG_41_left.data(:,:,i);
    data_to_warp_41_left(:,evLatency_left(i,end)+1:end) = [];
    warped_data_41_left{i} = warpmats_left{i} * data_to_warp_41_left';
end

% create new set
numEpochs = length(warped_data_39_left);
matrixSize = size(warped_data_39_left{1}');
inverted_data_39_left = zeros([matrixSize, numEpochs], 'single');
inverted_data_41_left = zeros([matrixSize, numEpochs], 'single');

for i = 1:numEpochs
    inverted_data_39_left(:, :, i) = warped_data_39_left{i}';
    inverted_data_41_left(:, :, i) = warped_data_41_left{i}';
end
EEG_39_warped_left = EEG_39_left;
EEG_39_warped_left.data = inverted_data_39_left;
EEG_41_warped_left = EEG_41_left;
EEG_41_warped_left.data = inverted_data_41_left;

% match latencies in set
for i = 1:length(EEG_39_left.epoch)
    for j = 2:length(newLatency_ms)+1
        EEG_39_warped_left.epoch(i).eventlatency{j} = [newLatency_ms(j-1)];
        EEG_41_warped_left.epoch(i).eventlatency{j} = [newLatency_ms(j-1)];
    end
end

% warp right turns
warped_data_39_right = cell(length(EEG_39_right.epoch),1);     % initialize cell array
warped_data_41_right = cell(length(EEG_41_right.epoch),1);
for i = 1:length(EEG_39_right.epoch)
    data_to_warp_39_right = EEG_39_right.data(:,:,i);
    data_to_warp_39_right(:,evLatency_right(i,end)+1:end) = [];
    warped_data_39_right{i} = warpmats_right{i} * data_to_warp_39_right';
    data_to_warp_41_right = EEG_41_right.data(:,:,i);
    data_to_warp_41_right(:,evLatency_right(i,end)+1:end) = [];
    warped_data_41_right{i} = warpmats_right{i} * data_to_warp_41_right';
end

% create new set
numEpochs = length(warped_data_39_right);
matrixSize = size(warped_data_39_right{1}');
inverted_data_39_right = zeros([matrixSize, numEpochs], 'single');
inverted_data_41_right = zeros([matrixSize, numEpochs], 'single');

for i = 1:numEpochs
    inverted_data_39_right(:, :, i) = warped_data_39_right{i}';
    inverted_data_41_right(:, :, i) = warped_data_41_right{i}';
end
EEG_39_warped_right = EEG_39_right;
EEG_39_warped_right.data = inverted_data_39_right;
EEG_41_warped_right = EEG_41_right;
EEG_41_warped_right.data = inverted_data_41_right;

% match latencies in set
for i = 1:length(EEG_39_right.epoch)
    for j = 2:length(newLatency_ms)+1
        EEG_39_warped_right.epoch(i).eventlatency{j} = [newLatency_ms(j-1)];
        EEG_41_warped_right.epoch(i).eventlatency{j} = [newLatency_ms(j-1)];
    end
end
fprintf('You have successfully timewarped the data!!\n');

%% save sets

% Define file paths and names
output_file_path = ['C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\3_EPOCHED\P' sub_id];

% Check if the folder exists, and create it if not
if ~exist(output_file_path, 'dir')
    mkdir(output_file_path);
end

% left
output_file_name = ['P' sub_id '_39' condition{iblock} '_left_epoched.set'];
pop_saveset(EEG_39_warped_left, 'filename', output_file_name, 'filepath', output_file_path);
fprintf('EEG data saved as %s\n', output_file_name);
output_file_name = ['P' sub_id '_41' condition{iblock} '_left_epoched.set'];
pop_saveset(EEG_41_warped_left, 'filename', output_file_name, 'filepath', output_file_path);
fprintf('EEG data saved as %s\n', output_file_name);

% right
output_file_name = ['P' sub_id '_39' condition{iblock} '_right_epoched.set'];
pop_saveset(EEG_39_warped_right, 'filename', output_file_name, 'filepath', output_file_path);
fprintf('EEG data saved as %s\n', output_file_name);
output_file_name = ['P' sub_id '_41' condition{iblock} '_right_epoched.set'];
pop_saveset(EEG_41_warped_right, 'filename', output_file_name, 'filepath', output_file_path);
fprintf('EEG data saved as %s\n', output_file_name);
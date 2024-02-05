%% calculate assrMI left

assrPOW39_left = cell(length(EEG_39_warped_left.epoch),1);      % Initialize array for 39 Hz power
assrPOW41_left = cell(length(EEG_41_warped_left.epoch),1);      % Initialize array for 41 Hz power

for itrial = 1:length(EEG_39_warped_left.epoch)
    assrPOW39_left{itrial} = abs(EEG_39_warped_left.data(:,:,itrial)).^2;
    assrPOW39_left{itrial} = mean(assrPOW39_left{itrial});
    assrPOW41_left{itrial} = abs(EEG_41_warped_left.data(:,:,itrial)).^2;
    assrPOW41_left{itrial} = mean(assrPOW41_left{itrial});
end

assrPOW39_left_double = zeros(length(EEG_39_warped_left.data(1,1,:)), length(EEG_39_warped_left.data(1,:,1)));
assrPOW41_left_double = zeros(length(EEG_41_warped_left.data(1,1,:)), length(EEG_41_warped_left.data(1,:,1)));
for i = 1:length(EEG_39_warped_left.epoch)
    assrPOW39_left_double(i,:) = assrPOW39_left{i};
    assrPOW41_left_double(i,:) = assrPOW41_left{i};
end

assrMI_left = zeros(length(assrPOW39_left_double(:,1)), length(assrPOW39_left_double(1,:)));
for j = 1:length(assrMI_left(:,1))
    for i = 1:length(assrPOW39_left_double)
        assrMI_left(j,i) = (assrPOW39_left_double(j,i) - assrPOW41_left_double(j,i)) / (abs(assrPOW39_left_double(j,i)) + abs(assrPOW41_left_double(j,i)));
    end
end

assrMI_left_mean = mean(assrMI_left,1);
assrMI_left_zscore = zscore(assrMI_left_mean);

%% calculate assrMI right

assrPOW39_right = cell(length(EEG_39_warped_right.epoch),1);      % Initialize array for 39 Hz power
assrPOW41_right = cell(length(EEG_41_warped_right.epoch),1);      % Initialize array for 41 Hz power

for itrial = 1:length(EEG_39_warped_right.epoch)
    assrPOW39_right{itrial} = abs(EEG_39_warped_right.data(:,:,itrial)).^2;
    assrPOW39_right{itrial} = mean(assrPOW39_right{itrial});
    assrPOW41_right{itrial} = abs(EEG_41_warped_right.data(:,:,itrial)).^2;
    assrPOW41_right{itrial} = mean(assrPOW41_right{itrial});
end

assrPOW39_right_double = zeros(length(EEG_39_warped_right.data(1,1,:)), length(EEG_39_warped_right.data(1,:,1)));
assrPOW41_right_double = zeros(length(EEG_41_warped_right.data(1,1,:)), length(EEG_41_warped_right.data(1,:,1)));
for i = 1:length(EEG_39_warped_right.epoch)
    assrPOW39_right_double(i,:) = assrPOW39_right{i};
    assrPOW41_right_double(i,:) = assrPOW41_right{i};
end

assrMI_right = zeros(length(assrPOW39_right_double(:,1)), length(assrPOW39_right_double(1,:)));
for j = 1:length(assrMI_right(:,1))
    for i = 1:length(assrPOW39_right_double)
        assrMI_right(j,i) = (assrPOW39_right_double(j,i) - assrPOW41_right_double(j,i)) / (abs(assrPOW39_right_double(j,i)) + abs(assrPOW41_right_double(j,i)));
    end
end

assrMI_right_mean = mean(assrMI_right,1);
assrMI_right_zscore = zscore(assrMI_right_mean);

%% save in data structures

% for frequencies individually for each direction

if ismember(subject, cell2mat(VR))
    % for ASSR Power
    if isub == 1
        VRsub = 1;
    end
    if isub == 2
        VRsub = 2;
    end
    if isub == 4
        VRsub = 3;
    end
    % data structure: data{sub, block}{left/right, 39/41}
    assrPOW_data_VR{VRsub, iblock}{1, 1} = assrPOW39_left;
    assrPOW_data_VR{VRsub, iblock}{1, 2} = assrPOW41_left;
    assrPOW_data_VR{VRsub, iblock}{2, 1} = assrPOW39_right;
    assrPOW_data_VR{VRsub, iblock}{2, 2} = assrPOW41_right;

    % for assrMI
    % data structure: data{sub, block}{left/right}
    assrMI_data_VR{VRsub, iblock}{1} = assrMI_left;
    assrMI_data_VR{VRsub, iblock}{2} = assrMI_right;

    % for z-scored assrMI
    % data structure: data{sub, block}{left/right}
    assrMI_data_zscore_VR{VRsub, iblock}{1} = assrMI_left_zscore;
    assrMI_data_zscore_VR{VRsub, iblock}{2} = assrMI_right_zscore;
else
    if isub == 3
        noVRsub = 1;
    end
    if isub == 5
        noVRsub = 2;
    end
    if isub == 6
        noVRsub = 3;
    end
    noVRindex = iblock - 2;
    % data structure: data{sub, block}{left/right, 39/41}
    assrPOW_data_noVR{noVRsub, noVRindex}{1, 1} = assrPOW39_left;
    assrPOW_data_noVR{noVRsub, noVRindex}{1, 2} = assrPOW41_left;
    assrPOW_data_noVR{noVRsub, noVRindex}{2, 1} = assrPOW39_right;
    assrPOW_data_noVR{noVRsub, noVRindex}{2, 2} = assrPOW41_right;

    % for assrMI
    % data structure: data{sub, block}{left/right}
    assrMI_data_noVR{noVRsub, noVRindex}{1} = assrMI_left;
    assrMI_data_noVR{noVRsub, noVRindex}{2} = assrMI_right;

    % for z-scored assrMI
    % data structure: data{sub, block}{left/right}
    assrMI_data_zscore_noVR{noVRsub, noVRindex}{1} = assrMI_left_zscore;
    assrMI_data_zscore_noVR{noVRsub, noVRindex}{2} = assrMI_right_zscore;
end


index_of_20percent = round(length(assrMI_left) * 0.2);
index_of_80percent = round(length(assrMI_left) * 0.8);

assrMI_left = assrMI_left(:, index_of_20percent:index_of_80percent);
assrMI_right = assrMI_right(:, index_of_20percent:index_of_80percent);
assrMI_left_mean = assrMI_left_mean(:, index_of_20percent:index_of_80percent);
assrMI_right_mean = assrMI_right_mean(:, index_of_20percent:index_of_80percent);
assrMI_left_zscore = assrMI_left_zscore(:, index_of_20percent:index_of_80percent);
assrMI_right_zscore = assrMI_right_zscore(:, index_of_20percent:index_of_80percent);
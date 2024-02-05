%% assrLI script

clear all;
close all;

%% what do you want to do?

figures                 = 'yes'; % if you want figure 'yes', if they already exist 'no'
statistics_noVR         = 'yes'; % if you want stas 'yes', if not 'no'
statistics_VR           = 'no'; % if you want stas 'yes', if not 'no'
statistics_pow          = 'no'; % if you want stas 'yes', if not 'no'

%% specifications

goodSubs                = {1, 2, 8, 9, 10, 11};
VR                      = {1, 2, 9};
noVR                    = {8, 10, 11};
condition               = {'_A_VR','_B_VR','_A','_B'};

subject = goodSubs{isub};
sub_id = sprintf('%03d', subject);

if ismember(subject, cell2mat(VR))
    blocks = 1:2;
else
    blocks = 3:4;
end

%% initiate storing cells

assrMI_left_avg_VR      = 0;
assrMI_left_avg_noVR    = 0;
assrMI_right_avg_VR     = 0;
assrMI_right_avg_noVR   = 0;
assrPOW_data_VR         = cell(length(VR), 2);
assrPOW_data_noVR       = cell(length(noVR), 2);
assrMI_data_VR          = cell(length(VR), 2);
assrMI_data_noVR        = cell(length(noVR), 2);
assrMI_data_zscore_VR   = cell(length(VR), 2);
assrMI_data_zscore_noVR = cell(length(noVR), 2);

%% start looping

for isub = 1:length(goodSubs)
    for iblock = blocks
        %% import data

        data_path       = ['C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\3_EPOCHED\P' sub_id];
        study_file_path = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\';

        output_file_name = ['P' sub_id '_39' condition{iblock} '_left_epoched.set'];
        EEG_39_warped_left = pop_loadset('filename', output_file_name, 'filepath', data_path);
        output_file_name = ['P' sub_id '_41' condition{iblock} '_left_epoched.set'];
        EEG_41_warped_left = pop_loadset('filename', output_file_name, 'filepath', data_path);

        output_file_name = ['P' sub_id '_39' condition{iblock} '_right_epoched.set'];
        EEG_39_warped_right = pop_loadset('filename', output_file_name, 'filepath', data_path);
        output_file_name = ['P' sub_id '_41' condition{iblock} '_right_epoched.set'];
        EEG_41_warped_right = pop_loadset('filename', output_file_name, 'filepath', data_path);

        fprintf('You have successfully imported the data!!\n');

        %% select lateralized frontal electrodes
        % this should be adapted to the needs of our specific data

        channels_to_keep = [5, 7];
        EEG_39_warped_left = pop_select(EEG_39_warped_left, 'channel', channels_to_keep);
        EEG_41_warped_left = pop_select(EEG_41_warped_left, 'channel', channels_to_keep);
        EEG_39_warped_right = pop_select(EEG_39_warped_right, 'channel', channels_to_keep);
        EEG_41_warped_right = pop_select(EEG_41_warped_right, 'channel', channels_to_keep);
        fprintf('You have successfully selected F3 and F4!!\n');

        %% ASSR LI calculation

        ASSR_LI_calc

        %% figures?

        if figures == 'yes'
            ASSR_LI_plot
        end

        %% statistics noVR?

        if statistics_noVR == 'yes'
            ASSR_LI_stat
        end

        %% statistics VR?

        if statistics_VR == 'yes'
            ASSR_LI_stat
        end

        %% separate power?

        if statistics_pow == 'yes'
            ASSR_LI_power
        end
    end
end



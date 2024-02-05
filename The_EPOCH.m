%% epoch and warp script 

clear all;
close all;

%% specifications

goodSubs                = {1, 2, 8, 9, 10, 11};
VR                      = {1, 2, 9};
noVR                    = {8, 10, 11};
condition               = {'_A_VR','_B_VR','_A','_B'};
study_file_path         = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\';
directions              = {'left', 'right'};
maxLatency_ms           = 0;

subject = goodSubs{isub};
sub_id = sprintf('%03d', subject);  
if ismember(subject, cell2mat(VR))
    blocks = 1:2;
else
    blocks = 3:4;
end

%% start looping
for isub = 1:length(goodSubs)       
    for iblock          = blocks
        %% import

        data_path       = [study_file_path '2_PREPROCESSED\P' sub_id];
        file_name_39    = ['P' sub_id condition{iblock} '_39_bandpass.set'];
        file_name_41    = ['P' sub_id condition{iblock} '_41_bandpass.set'];
        EEG_39          = pop_loadset('filename', file_name_39, 'filepath', data_path);
        EEG_41          = pop_loadset('filename', file_name_41, 'filepath', data_path);
        fprintf('You have successfully imported the data!!\n');

        %% find max durations

        ASSR_epoch_durations

        %% epoch data

        ASSR_epoch_exec

        %% set newLatency for timewarp

        ASSR_epoch_newLatency
    end
end

for isub = 1:length(goodSubs)       
    for iblock          = blocks
        %% import

        data_path       = [study_file_path '2_PREPROCESSED\P' sub_id];
        file_name_39    = ['P' sub_id condition{iblock} '_39_bandpass.set'];
        file_name_41    = ['P' sub_id condition{iblock} '_41_bandpass.set'];
        EEG_39          = pop_loadset('filename', file_name_39, 'filepath', data_path);
        EEG_41          = pop_loadset('filename', file_name_41, 'filepath', data_path);
        fprintf('You have successfully imported the data!!\n');

        %% find max durations

        ASSR_epoch_durations

        %% epoch data

        ASSR_epoch_exec

        %% set evLatency

        ASSR_epoch_evLatency
        
        %% warp it baby

        ASSR_epoch_warp
        
    end
end









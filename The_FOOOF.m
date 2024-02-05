%% FOOOF script 

clear all;
close all;

%% specifications

goodSubs                                = {1, 2, 8, 9, 10, 11};
VR                                      = {1, 2, 9};
noVR                                    = {8, 10, 11};
condition                               = {'_A_VR','_B_VR','_A','_B'};
study_file_path                         = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\';
direction                               = {'left', 'right'};

subject = goodSubs{isub};
sub_id = sprintf('%03d', subject);

if ismember(subject, cell2mat(VR))
    blocks = 1:2;
else
    blocks = 3:4;
end

%% initiate storing cells

fooof_results_VR                        = cell(length(VR), length(2));
fooof_results_noVR                      = cell(length(noVR), length(2));

%% start looping 

for isub = 1:length(goodSubs)       
    for iblock                          = blocks 
        %% import 

        data_path                       = [study_file_path '2_PREPROCESSED\P' sub_id];
        file_name                       = ['P' sub_id condition{iblock} '_preprocessed.set'];
        EEG                             = pop_loadset('filename', file_name, 'filepath', data_path);
        fprintf('You have successfully imported the data!!\n');

        %% select lateralized frontal electrodes

        channels_to_keep = [5, 7];
        EEG = pop_select(EEG, 'channel', channels_to_keep);
        fprintf('You have successfully selected F3 and F4!!\n');

        %% split sets

        data_left_turns                 = EEG.data;
        data_right_turns                = EEG.data;

        %% prepare data for FOOOF

        ASSR_fooof_prep

        %% execute FOOOF

        settings = [];
        settings.max_n_peaks            = 4;
        settings.peak_threshold         = 1;
        settings.peak_width_limts       = 1;
        % settings.max_n_peaks
        % settings.min_peak_height
        % settings.peak_threshold
        % settings.aperiodic_mode
        % settings.verbose

        ASSR_fooof_exec

        %% plot results and save plot

        ASSR_fooof_plot
    end
end

%% export table and leave matlab

ASSR_fooof_export

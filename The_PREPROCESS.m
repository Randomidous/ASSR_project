%% Preprocess script 

clear all;
close all;

%% specifications

condition                   = {'_A_VR','_B_VR','_A','_B'};
directions                  = {'left', 'right'};       
isub                        = 14;
sub_id                      = sprintf('%03d', isub);
iblock                      = 1;
study_file_path             = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\';
output_file_name            = ['P' sub_id condition{iblock} '.set'];
output_file_path            = [study_file_path '1_SETS\P' sub_id];

%% import data 

% Check if the .set file already exists
if exist(fullfile(output_file_path, output_file_name), 'file')
    % .set file already exists, load the existing set
    EEG = pop_loadset('filename', output_file_name, 'filepath', output_file_path);
    fprintf('Set successfully loaded!!\n');
else
    % .set file doesn't exist, perform XDF import and create a new set
    file_path = [study_file_path '0_RAW\P' sub_id '\Block' condition{iblock} '.xdf'];
    EEG = pop_loadxdf(file_path, 'streamname', 'EEGstream EE225', 'streamtype', 'EEG', 'exclude_markerstreams', {}, 'HandleJitterRemoval', 0);
    EEG.setname = ['P' sub_id condition{iblock}];

    % Save .set in 1_SETS
    if ~exist(output_file_path, 'dir')
        mkdir(output_file_path);
    end

    EEG = pop_chanedit(EEG);
    pop_saveset(EEG, 'filename', output_file_name, 'filepath', output_file_path);
    fprintf('EEG data saved as %s\n', output_file_name);
end

%% pre-process data 

% bandpass filter
EEG = pop_eegfiltnew(EEG, 'locutoff', 1,'hicutoff', 100, 'plotfreqz', 1);
% average reference
EEG = pop_reref(EEG, []);
% remove line noise
EEG = pop_zapline_plus(EEG, 'noisefreqs',[50 90],'coarseFreqDetectPowerDiff',4,'chunkLength',0,'adaptiveNremove',1,'fixedNremove',1);
% clean data 
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','BurstCriterion',20,'WindowCriterion','off', 'BurstRejection','off');
% run PCA and ICA
EEG = pop_runica(EEG, 'icatype', 'runica', 'pca', 40); % EEG = pop_runamica(EEG, 'maxiter', 2000, 'pcakeep', 40); % EEG = eeg_checkset(EEG); % EEG = pop_multifit(EEG, [], 'icatype', 'amica', 'dataset', 1, 'intra', 'no');
% label ICs
EEG = pop_iclabel(EEG, 'lite');

fprintf('You have successfully pre-processed the data!!\n');

%% inspect data 

% PSD per channel
load layout.mat
[Pxx,F] = pwelch(EEG.data',1000,[],1000,500,'power');
figure;
for ichan = 1:length(EEG.data(:,1))
    subplot(8,8,ichan)
    plot(F,Pxx(:,ichan))
    xlim([1 100])
    name = layout.label(ichan,1);
    title(name)
end

% per component
figure;
for icomp = 1:length(EEG.icaact(:,1))
    subplot(7,6,icomp)
    [Pxx,F] = pwelch(EEG.icaact(icomp,:),1000,[],1000,500,'power');
    plot(F,mean(Pxx,2));
    xlim([10 45])
    title(num2str(icomp))
end

% component topography
pop_topoplot(EEG, 0, 1:length(EEG.icaact(:,1)) ,['P' sub_id condition{iblock} '.set'],[6 7] ,0,'electrodes','on');

fprintf('Can you spot decent ASSRs?\n');

%% select ICs

user_input = input('Enter the components to keep (e.g., [12 14]): ');
if ~isempty(user_input)
    EEG = pop_subcomp(EEG, user_input, 0, 1);
    eeglab redraw
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname', output_file_name,'gui','off');
    fprintf('You have successfully discarded the specified ICs!\n');
else
    fprintf('No components selected. Continuing without discarding.\n');
end

%% inspect data

[Pxx,F] = pwelch(EEG.data',1000,[],1000,500,'power');
figure;
for ichan = 1:length(EEG.data(:,1))
    subplot(8,8,ichan)
    plot(F,Pxx(:,ichan))
    xlim([1,45])
    name = layout.label(ichan,1);
    title(name)
end

%% further steps

% remove events before and after exp block
% rename events
% combine events with same latency
% rename events again
% look for other combined events and delete them
% remove last heading event
% adjust latency so that first mid-point event = 0
% Check for consecutive "Turning" events and delete the second one

ASSR__pre_preparation.m

%% are you special? Extra attention needed?

if isub == 11
    ASSR_pre_extrawurst
end

%% Save the _preprocessed.set

output_file_path = [study_file_path '2_PREPROCESSED\P' sub_id];

% Check if the folder exists, and create it if not
if ~exist(output_file_path, 'dir')
    mkdir(output_file_path);
end

output_file_name = ['P' sub_id condition{iblock} '_preprocessed.set'];
pop_saveset(EEG, 'filename', output_file_name, 'filepath', output_file_path);
fprintf('EEG data saved as %s\n', output_file_name);

%% bandpassfilter around signals of interest

EEG_39 = EEG;
EEG_41 = EEG;
ban = 0.5;

EEG_39.data = ft_preproc_bandpassfilter(EEG.data, 500, [39-ban 39+ban], 4, 'but','twopass');
EEG_41.data = ft_preproc_bandpassfilter(EEG.data, 500, [41-ban 41+ban], 4, 'but','twopass');
fprintf('You have successfully applied bandpass filters!!\n');

%% hilbert transform

EEG_39_hilbert = EEG_39;
EEG_41_hilbert = EEG_41;

for i = 1:length(EEG_39.data(:,1))
    EEG_39_hilbert.data(i,:) = abs(hilbert(EEG_39.data(i,:)));
    EEG_41_hilbert.data(i,:) = abs(hilbert(EEG_41.data(i,:)));
end
fprintf('You have successfully hilbert-transformed the data!!\n');

%% Save the _preprocessed.set

% Define file paths and names
output_file_path = [study_file_path '2_PREPROCESSED\P' sub_id];

% Check if the folder exists, and create it if not
if ~exist(output_file_path, 'dir')
    mkdir(output_file_path);
end

output_file_name = ['P' sub_id condition{iblock} '_39_bandpass.set'];
pop_saveset(EEG_39_hilbert, 'filename', output_file_name, 'filepath', output_file_path);
fprintf('EEG data saved as %s\n', output_file_name);

output_file_name = ['P' sub_id condition{iblock} '_41_bandpass.set'];
pop_saveset(EEG_41_hilbert, 'filename', output_file_name, 'filepath', output_file_path);
fprintf('EEG data saved as %s\n', output_file_name);


%% now epoch

EEG_39_left = pop_epoch(EEG_39, {'Torso entered MidPoint:TurningLeft'}, [lower_window  upper_window], 'epochinfo', 'yes');
EEG_41_left = pop_epoch(EEG_41, {'Torso entered MidPoint:TurningLeft'}, [lower_window  upper_window], 'epochinfo', 'yes');
EEG_39_right = pop_epoch(EEG_39, {'Torso entered MidPoint:TurningRight'}, [lower_window  upper_window], 'epochinfo', 'yes');
EEG_41_right = pop_epoch(EEG_41, {'Torso entered MidPoint:TurningRight'}, [lower_window  upper_window], 'epochinfo', 'yes');

[~,~] = pop_autorej(EEG_39_left, 'eegplot', 'off', 'nogui', 'on');
[~,~] = pop_autorej(EEG_41_left, 'eegplot', 'off', 'nogui', 'on');
[~,~] = pop_autorej(EEG_39_right, 'eegplot', 'off', 'nogui', 'on');
[~,~] = pop_autorej(EEG_41_right, 'eegplot', 'off', 'nogui', 'on');

% Automatically identify and remove epochs with eventtype size smaller than 9

sketchy_epochs_left = find(arrayfun(@(x) numel(x.eventtype) < 8, EEG_39_left.epoch));
sketchy_epochs_right = find(arrayfun(@(x) numel(x.eventtype) < 8, EEG_39_right.epoch));

% remove identified epochs
EEG_39_left = pop_selectevent(EEG_39_left, 'omitepoch', sketchy_epochs_left ,'deleteevents','off','deleteepochs','on','invertepochs','off');
EEG_41_left = pop_selectevent(EEG_41_left, 'omitepoch', sketchy_epochs_left ,'deleteevents','off','deleteepochs','on','invertepochs','off');
EEG_39_right = pop_selectevent(EEG_39_right, 'omitepoch', sketchy_epochs_right ,'deleteevents','off','deleteepochs','on','invertepochs','off');
EEG_41_right = pop_selectevent(EEG_41_right, 'omitepoch', sketchy_epochs_right ,'deleteevents','off','deleteepochs','on','invertepochs','off');

% special attention for P011
if isub == 6 && iblock == 4
    evLatency_ms_left = zeros(length(EEG_39_left.epoch), length(EEG_39_left.epoch(1).eventlatency)); % create evLatency matrix
    evLatency_ms_right = zeros(length(EEG_39_right.epoch), length(EEG_39_right.epoch(1).eventlatency)); % create evLatency matrix


    for i = 1:length(EEG_39_left.epoch)
        for j = 1:length(EEG_39_left.epoch(i).eventlatency)
            evLatency_ms_left(i, j) = cell2mat(EEG_39_left.epoch(i).eventlatency(1,j));
        end
    end
    for i = 1:length(EEG_39_right.epoch)
        for j = 1:length(EEG_39_right.epoch(i).eventlatency)
            evLatency_ms_right(i, j) = cell2mat(EEG_39_right.epoch(i).eventlatency(1,j));
        end
    end

    sketchy_epochs_left = find(any(evLatency_ms_left < 0, 2));
    sketchy_epochs_right = find(any(evLatency_ms_right < 0, 2));


    % remove identified epochs
    if ~isempty(sketchy_epochs_left)
        EEG_39_left = pop_selectevent(EEG_39_left, 'omitepoch', sketchy_epochs_left ,'deleteevents','off','deleteepochs','on','invertepochs','off');
        EEG_41_left = pop_selectevent(EEG_41_left, 'omitepoch', sketchy_epochs_left ,'deleteevents','off','deleteepochs','on','invertepochs','off');
        EEG_39_right = pop_selectevent(EEG_39_right, 'omitepoch', sketchy_epochs_right ,'deleteevents','off','deleteepochs','on','invertepochs','off');
        EEG_41_right = pop_selectevent(EEG_41_right, 'omitepoch', sketchy_epochs_right ,'deleteevents','off','deleteepochs','on','invertepochs','off');

    else
        disp('No sketchy epochs found.');
    end
end
fprintf('You have successfully epoched the data!!\n');
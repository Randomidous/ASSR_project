%% set evLatency and newLatency

% create evLatency matrix (note: only necessary for once frequency, since they share events)
evLatency_ms_left = zeros(length(EEG_39_left.epoch), length(EEG_39_left.epoch(1).eventlatency)); % create evLatency matrix
evLatency_ms_right = zeros(length(EEG_39_right.epoch), length(EEG_39_right.epoch(1).eventlatency)); % create evLatency matrix

for i = 1:length(EEG_39_left.epoch)
    for j = 1:length(EEG_39_left.epoch(i).eventlatency)
        evLatency_ms_left(i, j) = cell2mat(EEG_39_left.epoch(i).eventlatency(1,j));
    end
end
evLatency_ms_left(:,10:end) = [];    % remove everything after trial end
evLatency_ms_left(:,1) = [];         % remove zero

for i = 1:length(EEG_39_right.epoch)
    for j = 1:length(EEG_39_right.epoch(i).eventlatency)
        evLatency_ms_right(i, j) = cell2mat(EEG_39_right.epoch(i).eventlatency(1,j));
    end
end
evLatency_ms_right(:,10:end) = [];    % remove everything after trial end
evLatency_ms_right(:,1) = [];         % remove zero
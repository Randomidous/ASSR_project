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

% if isub == 6 && iblock == 4
%     evLatency_ms_left(:,8) = [];
%     evLatency_ms_right(:,8) = [];
% end

% create newLatency vector
newLatency_ms_left = zeros(1, length(evLatency_ms_left(1,:)));
for i = 1:length(evLatency_ms_left(1,:))
    newLatency_ms_left(i) = median(evLatency_ms_left(:,i));
end

newLatency_ms_right = zeros(1, length(evLatency_ms_right(1,:)));
for i = 1:length(evLatency_ms_right(1,:))
    newLatency_ms_right(i) = median(evLatency_ms_right(:,i));
end

% pick the higher one for timewarping and update global

if max(newLatency_ms_left) > max(newLatency_ms_right)
    newLatency_ms = newLatency_ms_left;
    % Update global maximum
    if max(newLatency_ms) > max(maxLatency_ms)
        maxLatency_ms = newLatency_ms;
    end
else
    newLatency_ms = newLatency_ms_right;
    % Update global maximum
    if max(newLatency_ms) > max(maxLatency_ms)
        maxLatency_ms = newLatency_ms;
    end
end

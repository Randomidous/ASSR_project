%% instead of trimming, remove events before and after exp block

if contains(condition{iblock}, 'B')
    block = 'B';
else
    block = 'A';
end

% Find the indices and latencies
event_indices_start = find(strcmp({EEG.event.type}, ['trial:start;trialNo:1;heading:' block]));
matching_events = find(cellfun(@(x) contains(x, 'trial:end;trialNo:'), {EEG.event.type}));

if isub == 9 && block == 'B'
    matching_events = find(cellfun(@(x) contains(x, 'trial:end;trialNo:50'), {EEG.event.type}));
end

event_indices_end = matching_events(end);
start_latency = EEG.event(event_indices_start(1)).latency;
end_latency = EEG.event(event_indices_end(end)).latency;

events_to_remove = find([EEG.event.latency] < start_latency | [EEG.event.latency] > end_latency); % Find events outside the specified latency range
EEG.event(events_to_remove) = []; % Remove events from EEG.event
for i = 1:length(EEG.event) % Update event indices after removal
    EEG.event(i).urevent = i;
end
fprintf('You have successfully removed events outside the specified latency range!\n');

%% rename events

event_type_to_split = 'Torso entered MidPoint'; % define event
event_indices_start = find(strcmp({EEG.event.type}, event_type_to_split)); % Find indices of the events
new_midpoints = [];

% if Block B
if contains(condition{iblock}, 'B')
    % iterate through and rename; if even index = right, if uneven = left
    for i = 1:length(event_indices_start)
        original_midpoints = EEG.event(event_indices_start(i));
        heading_event = original_midpoints;
        if mod(i, 2) == 0
            heading_event.type = 'Torso entered MidPoint:TurningRight';
        else
            heading_event.type = 'Torso entered MidPoint:TurningLeft';
        end
        new_midpoints = [new_midpoints heading_event]; %#ok<*AGROW>
    end
    EEG.event(event_indices_start) = new_midpoints;
    fprintf('You have successfully renamed the events!!\n');
else
    % iterate through and rename; if even index = right, if uneven = left
    for i = 1:length(event_indices_start)
        original_midpoints = EEG.event(event_indices_start(i));
        heading_event = original_midpoints;
        if mod(i, 2) == 0
            heading_event.type = 'Torso entered MidPoint:TurningLeft';
        else
            heading_event.type = 'Torso entered MidPoint:TurningRight';
        end
        new_midpoints = [new_midpoints heading_event]; %#ok<*AGROW>
    end
    EEG.event(event_indices_start) = new_midpoints;
    fprintf('You have successfully renamed the events!!\n');
end

%% combine events with same latency

tolerance = 5;
unique_latencies = unique([EEG.event.latency]);
for i = 1:length(unique_latencies)
    indices = find(abs([EEG.event.latency] - unique_latencies(i)) <= tolerance);
    if length(indices) > 1
        combined_type = '';
        for j = 1:length(indices)
            combined_type = strcat(combined_type, EEG.event(indices(j)).type, ';');
        end
        combined_type = strcat(combined_type, 'combined');
        EEG.event(indices(1)).type = combined_type;
        EEG.event(indices(2:end)) = [];
    end
end
fprintf('You have successfully combined events with the same latencies!!\n');

%% rename events again

for i = 1:length(EEG.event)
    if contains(EEG.event(i).type, 'TurningRight')
        EEG.event(i).type = 'Torso entered MidPoint:TurningRight';
    elseif contains(EEG.event(i).type, 'TurningLeft')
        EEG.event(i).type = 'Torso entered MidPoint:TurningLeft';
    end
end
fprintf('You have successfully renamed the events again!!\n');

%% look for other combined events and delete them

for j = length(EEG.event):-1:1
    if contains(EEG.event(j).type, 'combined')
        EEG.event(j) = [];
    end
end
fprintf('Combined events have been deleted!\n');


%% remove last heading event (there is one generated the last time torso enters midpoint)

event_indices = find(startsWith({EEG.event.type}, 'Torso entered MidPoint'));
if ~isempty(event_indices)
    last_event_index = event_indices(end);
    EEG.event(last_event_index) = [];
    fprintf('You have successfully removed extra events!!\n');
else
    fprintf('No events matching the specified type were found.\n');
end

%% adjust latency so that first mid-point event = 0

if contains(condition{iblock}, 'B')
    event_indices = find(startsWith({EEG.event.type}, 'Torso entered MidPoint:TurningLeft')); % Get the first value in EEG.event.latency
    first_event_index = event_indices(1);
    first_latency = EEG.event(first_event_index).latency;

    for i = 1:length(EEG.event) % Loop through each event and adjust the latency
        EEG.event(i).latency = EEG.event(i).latency - first_latency;
    end
    fprintf('You have successfully set first mid-point latency to 0!!\n');
else
    event_indices = find(startsWith({EEG.event.type}, 'Torso entered MidPoint:TurningRight')); % Get the first value in EEG.event.latency
    first_event_index = event_indices(1);
    first_latency = EEG.event(first_event_index).latency;

    for i = 1:length(EEG.event) % Loop through each event and adjust the latency
        EEG.event(i).latency = EEG.event(i).latency - first_latency;
    end
    fprintf('You have successfully set first mid-point latency to 0!!\n');
end

%% Check for consecutive "Turning" events and delete the second one

for j = 2:length(EEG.event)
    if contains(EEG.event(j).type, 'Turning') && contains(EEG.event(j-1).type, 'Turning')
        % Check if consecutive events have the same type containing "Turning"
        % Delete the second event
        EEG.event(j) = [];

        % Rename events following the removed index
        for k = j:length(EEG.event)
            if contains(EEG.event(k).type, 'TurningRight')
                EEG.event(k).type = strrep(EEG.event(k).type, 'TurningRight', 'TurningLeft');
            elseif contains(EEG.event(k).type, 'TurningLeft')
                EEG.event(k).type = strrep(EEG.event(k).type, 'TurningLeft', 'TurningRight');
            end
        end
        break;
    end
end

fprintf('I am getting annoyed!!\n');
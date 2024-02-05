%% preparation for FOOOF

%% alright let's go

data_left_turns = EEG.data;
data_right_turns = EEG.data;

if contains(condition{iblock}, 'A') % this will be for blocks heading A
    % right turns = torso entered lower curve left until torso exited upper curve left
    % left turns = torso entered lower curve right until torso exited upper curve right

    event_types_left = {'Torso entered LowerCurveLeft', 'Torso exited UpperCurveLeft'};
    event_types_right = {'Torso entered LowerCurveRight', 'Torso exited UpperCurveRight'};
    event_indices_left = cell(1, length(event_types_left));
    event_indices_right = cell(1, length(event_types_right));
    event_latencies_left = cell(1, length(event_types_left));
    event_latencies_right = cell(1, length(event_types_right));

    for i = 1:length(event_types_left)
        % Process left events
        event_indices_left{i} = find(strcmp({EEG.event.type}, event_types_left{i}));
        event_latencies_left{i} = [EEG.event(event_indices_left{i}).latency];
        for j = 2:length(event_indices_left{i})
            if (event_indices_left{i}(j) - event_indices_left{i}(j-1)) > 25
                % if the difference is more than 25, add another latency
                new_index = round((event_latencies_left{i}(j) + event_latencies_left{i}(j-1)) / 2);
                event_latencies_left{i} = [event_latencies_left{i}(1:j-1), new_index, event_latencies_left{i}(j:end)];
            end
        end

        % Process right events
        event_indices_right{i} = find(strcmp({EEG.event.type}, event_types_right{i}));
        event_latencies_right{i} = [EEG.event(event_indices_right{i}).latency];
        for j = 2:length(event_indices_right{i})
            if (event_indices_right{i}(j) - event_indices_right{i}(j-1)) > 25
                % if the difference is more than 25, add another latency
                new_index = round((event_latencies_right{i}(j) + event_latencies_right{i}(j-1)) / 2);
                event_latencies_right{i} = [event_latencies_right{i}(1:j-1), new_index, event_latencies_right{i}(j:end)];
            end
        end
    end


    % Will do left turns first
    num_windows_left = length(event_latencies_left{1,2});
    data_left_turns_selected = cell(1, num_windows_left);

    for iwindows = 1:num_windows_left % iterate over each time window for left turns
        latency_start_left = event_latencies_left{1,1}(iwindows);
        latency_end_left = event_latencies_left{1,2}(iwindows);

        % Extract the EEG data for left turns within the specified time window
        data_left_turns_selected{iwindows} = data_left_turns(:, latency_start_left:latency_end_left);
    end
    concatenated_data_left_turns = cat(2, data_left_turns_selected{:});

    % and same for right
    num_windows_right = length(event_latencies_right{1,1});
    data_right_turns_selected = cell(1, num_windows_right);

    for iwindows = 1:num_windows_right % iterate over each time window for right turns
        latency_start_right = event_latencies_right{1,1}(iwindows);
        latency_end_right = event_latencies_right{1,2}(iwindows);

        % Extract the EEG data for right turns within the specified time window
        data_right_turns_selected{iwindows} = data_right_turns(:, latency_start_right:latency_end_right);
    end
    concatenated_data_right_turns = cat(2, data_right_turns_selected{:});

else
    % this is for blocks heading B

    event_types_left = {'Torso entered UpperCurveLeft', 'Torso exited LowerCurveLeft'};
    event_types_right = {'Torso entered UpperCurveRight', 'Torso exited LowerCurveRight'};
    event_indices_left = cell(1, length(event_types_left));
    event_indices_right = cell(1, length(event_types_right));
    event_latencies_left = cell(1, length(event_types_left));
    event_latencies_right = cell(1, length(event_types_right));

    for i = 1:length(event_types_left)
        event_indices_left{i} = find(strcmp({EEG.event.type}, event_types_left{i}));
        event_indices_right{i} = find(strcmp({EEG.event.type}, event_types_right{i}));
        event_latencies_left{i} = [EEG.event(event_indices_left{i}).latency];
        event_latencies_right{i} = [EEG.event(event_indices_right{i}).latency];
    end

    % Will do left turns first
    num_windows_left = length(event_latencies_left{1,2});
    data_left_turns_selected = cell(1, num_windows_left);

    for iwindows = 1:num_windows_left % iterate over each time window for left turns
        latency_start_left = event_latencies_left{1,1}(iwindows);
        latency_end_left = event_latencies_left{1,2}(iwindows);

        % Extract the EEG data for left turns within the specified time window
        data_left_turns_selected{iwindows} = data_left_turns(:, latency_start_left:latency_end_left);
    end
    concatenated_data_left_turns = cat(2, data_left_turns_selected{:});

    % and same for right
    num_windows_right = length(event_latencies_right{1,2});
    data_right_turns_selected = cell(1, num_windows_right);
    if isub == 6 && iblock == 4
        num_windows_right = num_windows_right-2;
    end

    for iwindows = 1:num_windows_right % iterate over each time window for right turns
        latency_start_right = event_latencies_right{1,1}(iwindows);
        latency_end_right = event_latencies_right{1,2}(iwindows);

        % Extract the EEG data for right turns within the specified time window
        data_right_turns_selected{iwindows} = data_right_turns(:, latency_start_right:latency_end_right);
    end
    concatenated_data_right_turns = cat(2, data_right_turns_selected{:});
end
%% find max durations

if contains(condition{iblock}, 'B')
    % Block B
    event_types = {'Torso entered MidPoint:TurningLeft', 'Torso entered MidPoint:TurningRight'};
    event_indices = cell(1, length(event_types));
    event_latencies = cell(1, length(event_types));

    for i = 1:length(event_types)
        event_indices{i} = find(strcmp({EEG_39.event.type}, event_types{i}));
        event_latencies{i} = [EEG_39.event(event_indices{i}).latency];
    end

    TurningLeft = event_latencies{1};
    TurningRight = event_latencies{2};

    % Remove the last entry from the array with more entries
    num_turning_right = length(TurningRight);
    num_turning_left = length(TurningLeft);
    if num_turning_right > num_turning_left
        diff = num_turning_right - num_turning_left;
        TurningRight(end-diff+1:end) = [];
    elseif num_turning_left > num_turning_right
        diff = num_turning_left - num_turning_right;
        TurningLeft(end-diff+1:end) = [];
    end

    durations_left =  TurningRight - TurningLeft;
    TurningLeft(1) = [];
    TurningRight(end) = [];
    durations_right = TurningLeft - TurningRight;

    Q1 = quantile(durations_left, 0.25);     % calculate inter-quartile range
    Q3 = quantile(durations_left, 0.75);
    IQR = Q3 - Q1;
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;
    outlier_indices_left = find(durations_left < lower_bound | durations_left > upper_bound);
    durations_without_outliers_left = durations_left;             % Remove outliers from the data
    durations_without_outliers_left(outlier_indices_left) = [];

    Q1 = quantile(durations_right, 0.25);     % calculate inter-quartile range
    Q3 = quantile(durations_right, 0.75);
    IQR = Q3 - Q1;
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;
    outlier_indices_right = find(durations_right < lower_bound | durations_right > upper_bound);
    durations_without_outliers_right = durations_right;             % Remove outliers from the data
    durations_without_outliers_right(outlier_indices_right) = [];

    % calculate max duration
    max_duration_left = max(durations_without_outliers_left);
    max_duration_sec_left = max_duration_left/500;
    fprintf('The max duration of a left turning trial is %.2f samples\n', max_duration_left);
    fprintf('That is %.2f seconds\n', max_duration_sec_left);

    max_duration_right = max(durations_without_outliers_right);
    max_duration_sec_right = max_duration_right/500;
    fprintf('The max duration of right turning a trial is %.2f samples\n', max_duration_right);
    fprintf('That is %.2f seconds\n', max_duration_sec_right);

    if max_duration_sec_left > max_duration_sec_right
        max_duration_sec = max_duration_sec_left;
        max_duration = max_duration_left;
    else
        max_duration_sec = max_duration_sec_right;
        max_duration = max_duration_right;
    end

    % set upper and lower window for epoching
    lower_window = -1;
    upper_window = max_duration_sec +1;

else
    % Block A
    event_types = {'Torso entered MidPoint:TurningRight', 'Torso entered MidPoint:TurningLeft'};
    event_indices = cell(1, length(event_types));
    event_latencies = cell(1, length(event_types));

    for i = 1:length(event_types)
        event_indices{i} = find(strcmp({EEG_39.event.type}, event_types{i}));
        event_latencies{i} = [EEG_39.event(event_indices{i}).latency];
    end

    TurningRight = event_latencies{1};
    TurningLeft = event_latencies{2};

    % Remove the last entry from the array with more entries
    num_turning_right = length(TurningRight);
    num_turning_left = length(TurningLeft);
    if num_turning_right > num_turning_left
        diff = num_turning_right - num_turning_left;
        TurningRight(end-diff+1:end) = [];
    elseif num_turning_left > num_turning_right
        diff = num_turning_left - num_turning_right;
        TurningLeft(end-diff+1:end) = [];
    end

    durations_right =  TurningLeft - TurningRight;
    % if isub == 8
    %     TurningRight = TurningRight(5:end)-13;
    % end
    TurningRight(1) = [];
    TurningLeft(end) = [];
    durations_left = TurningRight - TurningLeft;

    Q1 = quantile(durations_left, 0.25);     % calculate inter-quartile range
    Q3 = quantile(durations_left, 0.75);
    IQR = Q3 - Q1;
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;
    outlier_indices_left = find(durations_left < lower_bound | durations_left > upper_bound);
    durations_without_outliers_left = durations_left;             % Remove outliers from the data
    durations_without_outliers_left(outlier_indices_left) = [];

    Q1 = quantile(durations_right, 0.25);     % calculate inter-quartile range
    Q3 = quantile(durations_right, 0.75);
    IQR = Q3 - Q1;
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;
    outlier_indices_right = find(durations_right < lower_bound | durations_right > upper_bound);
    durations_without_outliers_right = durations_right;             % Remove outliers from the data
    durations_without_outliers_right(outlier_indices_right) = [];

    % calculate max duration
    max_duration_left = max(durations_without_outliers_left);
    max_duration_sec_left = max_duration_left/500;
    fprintf('The max duration of a left turning trial is %.2f samples\n', max_duration_left);
    fprintf('That is %.2f seconds\n', max_duration_sec_left);

    max_duration_right = max(durations_without_outliers_right);
    max_duration_sec_right = max_duration_right/500;
    fprintf('The max duration of a right turning trial is %.2f samples\n', max_duration_right);
    fprintf('That is %.2f seconds\n', max_duration_sec_right);

    if max_duration_sec_left > max_duration_sec_right
        max_duration_sec = max_duration_sec_left;
        max_duration = max_duration_left;
    else
        max_duration_sec = max_duration_sec_right;
        max_duration = max_duration_right;
    end

    % set upper and lower window for epoching
    lower_window = -1;
    upper_window = max_duration_sec +1;
end
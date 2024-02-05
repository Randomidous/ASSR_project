%% are you subject 11?

% Check for consecutive "Turning" events and delete the second one
    for often = 1:5
        j = 2;
        while j <= length(EEG.event)
            if contains(EEG.event(j).type, 'Turning') && contains(EEG.event(j-1).type, 'Turning')
                % Delete the second event
                EEG.event(j) = [];
                % Adjust j to stay at the same logical position after deletion
                continue; % Skip the increment of j to recheck at this position
            end
            j = j + 1;
        end

        % Check for "Turning" events 2 apart and delete the second one
        j = 3;
        while j <= length(EEG.event)
            if contains(EEG.event(j).type, 'Turning') && contains(EEG.event(j-2).type, 'Turning')
                % Delete the second event and the one immediately following it, if possible
                EEG.event(j) = [];
                if j <= length(EEG.event) % Ensure we do not exceed bounds
                    EEG.event(j) = []; % Adjust for potentially shifted indices
                end
                continue; % Skip the increment of j to recheck at this position
            end
            j = j + 1;
        end

        % Check for consecutive "Torso exited MidPoint" events and delete the latest one
        j = length(EEG.event);
        while j > 1
            if contains(EEG.event(j).type, 'Torso exited MidPoint') && contains(EEG.event(j-1).type, 'Torso exited MidPoint')
                EEG.event(j) = [];
                % Once a deletion is made, break since only one pair is targeted per iteration
                break;
            end
            j = j - 1;
        end
    end

    fprintf('You have successfully accounted for doubled mid-point events!!\n');
end

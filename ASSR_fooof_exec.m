%% this is where the FOOOF happens

fs = EEG.srate;         % sampling frequency in Hz
window_length = fs;     % 1 second window
overlap = fs/2;         % 0.5 second overlap
f_range = [1 45];      % change from Chen et al.'s range [1 48]

if ismember(subject, cell2mat(VR))
    if subject == 1
        VRsub = 1;
    end
    if subject == 2
        VRsub = 2;
    end
    if subject == 9
        VRsub = 3;
    end
    % left turns first
    [Pxx5, ~] = pwelch(concatenated_data_left_turns(1,:), hamming(window_length), overlap, [], fs, 'power'); % I have no idea why I need to transpose the data
    [Pxx7, F] = pwelch(concatenated_data_left_turns(2,:), hamming(window_length), overlap, [], fs, 'power'); % I have no idea why I need to transpose the data
    Pxx = (Pxx5 + Pxx7)./2; % average over two lateralized frontal electrodes
    fooof_results_VR{VRsub, iblock}{1,1} = fooof(F, Pxx, f_range, settings, true);

    % and le right turns
    [Pxx5, ~] = pwelch(concatenated_data_right_turns(1,:), hamming(window_length), overlap, [], fs, 'power'); % I have no idea why I need to transpose the data
    [Pxx7, F] = pwelch(concatenated_data_right_turns(2,:), hamming(window_length), overlap, [], fs, 'power'); % I have no idea why I need to transpose the data
    Pxx = (Pxx5 + Pxx7)./2; % average over two lateralized frontal electrodes
    fooof_results_VR{VRsub, iblock}{1,2} = fooof(F, Pxx, f_range, settings, true);
else
    % no VR
    if subject == 8
        noVRsub = 1;
    end
    if subject == 10
        noVRsub = 2;
    end
    if subject == 11
        noVRsub = 3;
    end
    noVRindex = iblock - 2;
    % left turns first
    [Pxx5, ~] = pwelch(concatenated_data_left_turns(1,:), hamming(window_length), overlap, [], fs, 'power'); % I have no idea why I need to transpose the data
    [Pxx7, F] = pwelch(concatenated_data_left_turns(2,:), hamming(window_length), overlap, [], fs, 'power'); % I have no idea why I need to transpose the data
    Pxx = (Pxx5 + Pxx7)./2; % average over two lateralized frontal electrodes
    fooof_results_noVR{noVRsub, noVRindex}{1,1} = fooof(F, Pxx, f_range, settings, true);

    % and le right turns
    [Pxx5, ~] = pwelch(concatenated_data_right_turns(1,:), hamming(window_length), overlap, [], fs, 'power'); % I have no idea why I need to transpose the data
    [Pxx7, F] = pwelch(concatenated_data_right_turns(2,:), hamming(window_length), overlap, [], fs, 'power'); % I have no idea why I need to transpose the data
    Pxx = (Pxx5 + Pxx7)./2; % average over two lateralized frontal electrodes
    fooof_results_noVR{noVRsub, noVRindex}{1,2} = fooof(F, Pxx, f_range, settings, true);
end

fprintf('I am become FOOOF. Cleaner of freqs.\n');
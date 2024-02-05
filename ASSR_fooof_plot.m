%% plot plot plot

% let's cut the crap
if ismember(subject, cell2mat(VR))
    left = fooof_results_VR{VRsub, iblock}{1, 1};
    right = fooof_results_VR{VRsub, iblock}{1, 2};
else
    left = fooof_results_noVR{noVRsub, noVRindex}{1, 1};
    right = fooof_results_noVR{noVRsub, noVRindex}{1, 2};
end

hFig = figure; % Create a figure
screenSize = get(0, 'ScreenSize'); % Get the screen size
set(hFig, 'Position', screenSize);

% PSD before and after FOOOF
subplot(3, 4, 1)
[Pxx, F] = pwelch(concatenated_data_left_turns(1, :), 1000, [], 1000, 500, 'power');
plot(F, mean(Pxx, 2));
xlim([30 45])
title('Left Turns - Using pwelsh');

subplot(3, 4, 2)
[Pxx, F] = pwelch(concatenated_data_right_turns(1, :), 1000, [], 1000, 500, 'power');
plot(F, mean(Pxx, 2));
xlim([30 45])
title('Right Turns - Using pwelsh');

subplot(3, 4, 3)
plot(left.freqs, left.power_spectrum);
xlim([30 45])
title('Left Turns - Using FOOOF');

subplot(3, 4, 4)
plot(right.freqs, right.power_spectrum);
xlim([30 45])
title('Right Turns - Using FOOOF');

% Aperiodic Fit and Oscillatory Peak
subplot(3, 4, 5)
plot(left.freqs, left.ap_fit, 'g', 'LineWidth', 2);
hold on;
scatter(left.peak_params(:, 1), left.peak_params(:, 2), 'ro');
xlabel('Frequency (Hz)');
ylabel('Power');
legend('Aperiodic Fit', 'Oscillatory Peaks', 'Location', 'southwest');
title('Aperiodic Fit and Oscillatory Peaks');

subplot(3, 4, 6)
plot(right.freqs, right.ap_fit, 'g', 'LineWidth', 2);
hold on;
scatter(right.peak_params(:, 1), right.peak_params(:, 2), 'ro');
xlabel('Frequency (Hz)');
ylabel('Power');
legend('Aperiodic Fit', 'Oscillatory Peaks', 'Location', 'southwest');
title('Aperiodic Fit and Oscillatory Peaks');

% Power Spectrum and FOOOF Model Fit
subplot(3, 4, 7)
plot(left.freqs, left.power_spectrum, 'b', 'LineWidth', 2);
hold on;
plot(left.freqs, left.fooofed_spectrum, 'r', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power');
legend('Original Spectrum', 'FOOOF Model Fit', 'Location', 'southwest');
title('Left Turns Power Spectrum and FOOOF Model Fit');

subplot(3, 4, 8)
plot(right.freqs, right.power_spectrum, 'b', 'LineWidth', 2);
hold on;
plot(right.freqs, right.fooofed_spectrum, 'r', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power');
legend('Original Spectrum', 'FOOOF Model Fit', 'Location', 'southwest');
title('Right Turns Power Spectrum and FOOOF Model Fit');

tightfig;

% and save
output_file_path = ['C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\P' sub_id];
if ~exist(output_file_path, 'dir')      % Check if the folder exists, and create it if not
    mkdir(output_file_path);
end

savename = ['C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\P' sub_id '\P' sub_id condition{iblock} '_fooof.png'];
saveas(hFig, savename);

close all
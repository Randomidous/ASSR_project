%% plot assrMI

hFig = figure; % Create a figure
screenSize = get(0, 'ScreenSize'); % Get the screen size
set(hFig, 'Position', screenSize);

% assrMI left
subplot(3, 2, 1)
x_values = 1:length(assrMI_left);
y_values = assrMI_left;
plot(x_values, y_values,  'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
xlim([0 length(assrMI_left)])
ylim([-1 1])
xlabel('sample');
ylabel('Value');
title(['P' sub_id ': assrMI turning left block ' condition{iblock}]);
grid on;
hold on
x_values = 1:length(assrMI_left);
y_values = assrMI_left_mean;
plot(x_values, y_values, 'LineWidth', 2, 'Color', 'red');
hold off

% assrMI right
subplot(3, 2, 2)
x_values = 1:length(assrMI_right);
y_values = assrMI_right;
plot(x_values, y_values,  'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
xlim([0 length(assrMI_right)])
ylim([-1 1])
xlabel('sample');
ylabel('Value');
title(['P' sub_id ': assrMI turning right block ' condition{iblock}]);
grid on;
hold on
x_values = 1:length(assrMI_right);
y_values = assrMI_right_mean;
plot(x_values, y_values, 'LineWidth', 2, 'Color', 'red');
hold off

% assrMI left zscore
subplot(3, 2, 3)
x_values = 1:length(assrMI_left_zscore);
y_values = assrMI_left_zscore;
plot(x_values, y_values,  'LineWidth', 2, 'Color', 'red');
xlim([0 length(assrMI_left_zscore)])
xlabel('sample');
ylabel('Value');
title(['P' sub_id ': z-scored assrMI turning left block ' condition{iblock}]);
grid on;

% assrMI right zscore
subplot(3, 2, 4)
x_values = 1:length(assrMI_right_zscore);
y_values = assrMI_right_zscore;
plot(x_values, y_values,  'LineWidth', 2, 'Color', 'red');
xlim([0 length(assrMI_right_zscore)])
xlabel('sample');
ylabel('Value');
title(['P' sub_id ': z-scored assrMI turning right block ' condition{iblock}]);
grid on;

% check normal distribution
subplot(3, 2, 5)
histogram(assrMI_left, 'Normalization', 'probability', 'EdgeColor', 'w');
title(['P' sub_id ': Histogram turning left block ' condition{iblock}]);

subplot(3, 2, 6)
histogram(assrMI_right, 'Normalization', 'probability', 'EdgeColor', 'w');
title(['P' sub_id ': Histogram turning right block ' condition{iblock}]);

% tightfig;

% and save
output_file_path = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\new';
if ~exist(output_file_path, 'dir')      % Check if the folder exists, and create it if not
    mkdir(output_file_path);
end

savename = ['C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\new\P' sub_id condition{iblock} '_assrLI_new.png'];
saveas(hFig, savename);

%% Prepare plot for seperate powers

% zscore
assrPOW39_left_zscore = zscore(mean(assrPOW39_left_double, 1));
assrPOW39_right_zscore = zscore(mean(assrPOW39_right_double, 1));
assrPOW41_left_zscore = zscore(mean(assrPOW41_left_double, 1));
assrPOW41_right_zscore = zscore(mean(assrPOW41_right_double, 1));

assrPOW39_left_zscore = assrPOW39_left_zscore(:, index_of_20percent:index_of_80percent);
assrPOW39_right_zscore = assrPOW39_right_zscore(:, index_of_20percent:index_of_80percent);
assrPOW41_left_zscore = assrPOW41_left_zscore(:, index_of_20percent:index_of_80percent);
assrPOW41_right_zscore = assrPOW41_right_zscore(:, index_of_20percent:index_of_80percent);

% Find the minimum value for each variable
min_values_left = min([assrPOW39_left_zscore, assrPOW41_left_zscore]);
min_values_right = min([assrPOW39_right_zscore, assrPOW41_right_zscore]);

% Subtract the minimum value from each variable
assrPOW39_left_zscore = assrPOW39_left_zscore - min(min_values_left);
assrPOW39_right_zscore = assrPOW39_right_zscore - min(min_values_right);
assrPOW41_left_zscore = assrPOW41_left_zscore - min(min_values_left);
assrPOW41_right_zscore = assrPOW41_right_zscore - min(min_values_right);

% ASSR Powers left
subplot(2, 2, 1)
x_values_left = 1:length(mean(assrPOW39_left_double(:, index_of_20percent:index_of_80percent), 1));
y_values_left = assrPOW39_left_double(:, index_of_20percent:index_of_80percent);
plot(x_values_left, y_values_left, 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
xlim([0 length(mean(assrPOW39_left_double(:, index_of_20percent:index_of_80percent), 1))])
ylim([-0.5 0.5])
hold on;
x_values_left = 1:length(mean(assrPOW39_left_double(:, index_of_20percent:index_of_80percent), 1));
y_values_left = mean(assrPOW39_left_double(:, index_of_20percent:index_of_80percent), 1);
plot(x_values_left, y_values_left, 'LineWidth', 2, 'Color', 'red');
hold on;
x_values_left = 1:length(mean(assrPOW39_left_double(:, index_of_20percent:index_of_80percent), 1));
y_values_left = -assrPOW41_left_double(:, index_of_20percent:index_of_80percent);
plot(x_values_left, y_values_left, 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
hold on;
x_values_left = 1:length(mean(assrPOW41_left_double(:, index_of_20percent:index_of_80percent), 1));
y_values_left = -mean(assrPOW41_left_double(:, index_of_20percent:index_of_80percent), 1);
plot(x_values_left, y_values_left, 'LineWidth', 2, 'Color', 'red');
grid on;
title(['P' sub_id ': ASSR powers turning left block ' condition{iblock}]);
hold off;

% ASSR Powers right
subplot(2, 2, 2)
x_values_right = 1:length(mean(assrPOW39_right_double(:, index_of_20percent:index_of_80percent), 1));
y_values_right = assrPOW39_right_double(:, index_of_20percent:index_of_80percent);
plot(x_values_right, y_values_right, 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
xlim([0 length(mean(assrPOW39_right_double(:, index_of_20percent:index_of_80percent), 1))])
ylim([-0.5 0.5])
hold on;
x_values_right = 1:length(mean(assrPOW39_right_double(:, index_of_20percent:index_of_80percent), 1));
y_values_right = mean(assrPOW39_right_double(:, index_of_20percent:index_of_80percent), 1);
plot(x_values_right, y_values_right, 'LineWidth', 2, 'Color', 'red');
hold on;
x_values_right = 1:length(mean(assrPOW39_right_double(:, index_of_20percent:index_of_80percent), 1));
y_values_right = -assrPOW41_right_double(:, index_of_20percent:index_of_80percent);
plot(x_values_right, y_values_right, 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
hold on;
x_values_right = 1:length(mean(assrPOW41_right_double(:, index_of_20percent:index_of_80percent), 1));
y_values_right = -mean(assrPOW41_right_double(:, index_of_20percent:index_of_80percent), 1);
plot(x_values_right, y_values_right, 'LineWidth', 2, 'Color', 'red');
grid on;
title(['P' sub_id ': ASSR powers turning right block ' condition{iblock}]);
hold off;

% ASSR Powers left zscore
subplot(2, 2, 3)
x_values_left = 1:length(mean(assrPOW39_left_zscore, 1));
y_values_left = assrPOW39_left_zscore;
plot(x_values_left, y_values_left, 'LineWidth', 2, 'Color', 'red');
xlim([0 length(mean(assrPOW39_left_zscore, 1))])
hold on;
x_values_left = 1:length(mean(assrPOW39_left_zscore, 1));
y_values_left = -assrPOW41_left_zscore;
plot(x_values_left, y_values_left, 'LineWidth', 2, 'Color', 'red');
hold on;
grid on;
title(['P' sub_id ': ASSR powers turning left block ' condition{iblock}]);
hold off;

% ASSR Powers right zscore
subplot(2, 2, 4)
x_values_right = 1:length(mean(assrPOW39_right_zscore, 1));
y_values_right = assrPOW39_right_zscore;
plot(x_values_right, y_values_right, 'LineWidth', 2, 'Color', 'red');
xlim([0 length(mean(assrPOW39_right_zscore, 1))])
hold on;
x_values_right = 1:length(mean(assrPOW39_right_zscore, 1));
y_values_right = -assrPOW41_right_zscore;
plot(x_values_right, y_values_right, 'LineWidth', 2, 'Color', 'red');
hold on;
grid on;
title(['P' sub_id ': ASSR powers turning right block ' condition{iblock}]);
hold off;

% tightfig;

% and save
output_file_path = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\new';
if ~exist(output_file_path, 'dir')      % Check if the folder exists, and create it if not
    mkdir(output_file_path);
end

savename = ['C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\new\P' sub_id condition{iblock} '_assrPOW_new.png'];
saveas(hFig, savename);

close all;
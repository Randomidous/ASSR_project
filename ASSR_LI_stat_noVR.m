%% get data

left_assrmi_VR = [];
right_assrmi_VR = [];
left_assrmi_noVR = [];
right_assrmi_noVR = [];
left_assrmi_noVR_meany = [];
right_assrmi_noVR_meany = [];

% data structure: data{sub, block}{left/right}
for i = 1:length(assrMI_data_VR(:,1))
    for j = 1:length(assrMI_data_VR(1,:))
        left_assrmi_VR = [left_assrmi_VR; assrMI_data_VR{i,j}{1,1}];
        right_assrmi_VR = [right_assrmi_VR; assrMI_data_VR{i,j}{1,2}];
        left_assrmi_noVR = [left_assrmi_noVR; assrMI_data_noVR{i,j}{1,1}];
        right_assrmi_noVR = [right_assrmi_noVR; assrMI_data_noVR{i,j}{1,2}];
        left_assrmi_noVR_meany = [left_assrmi_noVR_meany; mean(assrMI_data_noVR{i,j}{1,1},1)];
        right_assrmi_noVR_meany = [right_assrmi_noVR_meany; mean(assrMI_data_noVR{i,j}{1,2},1)];
    end
end

%% take middle 60% to only look at curves 

index_of_20percent = round(length(left_assrmi_VR) * 0.20); % same indices can be used for all data, since length should be the same after timewarping
index_of_80percent = round(length(left_assrmi_VR) * 0.80); % same indices can be used for all data, since length should be the same after timewarping

left_assrmi_noVR = left_assrmi_noVR(:,index_of_20percent:index_of_80percent);
right_assrmi_noVR = right_assrmi_noVR(:,index_of_20percent:index_of_80percent);
left_assrmi_noVR_meany = left_assrmi_noVR_meany(:,index_of_20percent:index_of_80percent);
right_assrmi_noVR_meany = right_assrmi_noVR_meany(:,index_of_20percent:index_of_80percent);

left_assrmi_noVR_mean = mean(left_assrmi_noVR, 1);
right_assrmi_noVR_mean = mean(right_assrmi_noVR, 1);
left_assrmi_noVR_zscore = zscore(left_assrmi_noVR_mean);
right_assrmi_noVR_zscore = zscore(right_assrmi_noVR_mean);


%% statistics: t-test for no VR

left_assrmi_noVR_mean2              = mean(left_assrmi_noVR, 2);
right_assrmi_noVR_mean2             = mean(right_assrmi_noVR, 2);
right_assrmi_noVR_mean2_fitting     = right_assrmi_noVR_mean2(1:length(left_assrmi_noVR_mean2),:); % match size
left_mean_mean                      = mean(left_assrmi_noVR_mean2);
right_mean_mean                     = mean(right_assrmi_noVR_mean2);

left                                = left_assrmi_noVR_mean2;
right                               = right_assrmi_noVR_mean2_fitting;

left = mean(left_assrmi_noVR_meany,2);
right = mean(right_assrmi_noVR_meany,2);

[h, p, ci, stats] = ttest(left, right, 'tail', 'right', 'alpha', 0.05);
ttest_noVR.h = h; 
ttest_noVR.p = p;
ttest_noVR.ci = ci;
ttest_noVR.stats = stats;

if p > 0.05
    fprintf('The difference between left and right assrMI values did not reach statistical significance (p = %f)!!\n', ttest_noVR.p);
else
    fprintf('The difference between left and right assrMI values did reach statistical significance (p = %f)!!\n', ttest_noVR.p);
end

mean_left_ttest = mean(left_assrmi_noVR_mean2);
std_left_ttest = std(left_assrmi_noVR_mean2);
mean_right_ttest = mean(right_assrmi_noVR_mean2_fitting);
std_right_ttest = std(right_assrmi_noVR_mean2_fitting);

%% permutation

% create struct for fieldtrip
datastr                     = [];
datastr.individual          = zeros(length(left_assrmi_noVR_meany(:,1)), 1, length(left_assrmi_noVR_meany));
datastr.time                = EEG_39_warped_left.times(:,index_of_20percent:index_of_80percent);
datastr.label               = {'N1'};
% datastr.dimord            = 'subj_chan_time';
left                        = datastr;
left.individual(:,1,:)      = left_assrmi_noVR_meany;
right                       = datastr;
right.individual(:,1,:)     = right_assrmi_noVR_meany;

% run test
lengths                     = length(left_assrmi_noVR_meany(:,1));
cfg                         = [];
cfg.parameter               = 'individual';
cfg.latency                 = 'all'; %[-2 2];
cfg.method                  = 'montecarlo';
cfg.statistic               = 'depsamplesT';
cfg.alpha                   = 0.05;
cfg.correcttail             = 'prob';
cfg.numrandomization        = 1000;
cfg.correctm                = 'cluster';
cfg.design(1:2*lengths, 1)  = [ones(lengths, 1); 2*ones(lengths, 1)];
cfg.design(1:2*lengths, 2)  = [repelem(1:6, 1, 1)'; repelem(1:6, 1, 1)'];
cfg.design                  = cfg.design';
cfg.ivar                    = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                    = 2; % the 2nd row in cfg.design contains the subject number
stat_noVRcorrect            = ft_timelockstatistics(cfg, left, right);

%% get non-1 areas

% Assuming stat_noVR.prob is your variable
nonOneIndices = find(stat_noVRcorrect.prob ~= 1);

% Alternatively, if you want to save it in a new variable
indicesNotEqualToOne = nonOneIndices;

lower_edge = min(indicesNotEqualToOne);
upper_edge = max(indicesNotEqualToOne);

%% plot assrMI grand mean

hFig = figure; % Create a figure
screenSize = get(0, 'ScreenSize'); % Get the screen size
set(hFig, 'Position', screenSize);

% assrMI left
subplot(2, 2, 1)
x_values = 1:length(left_assrmi_noVR_meany);
y_values = left_assrmi_noVR_meany;
plot(x_values, y_values,  'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
xlim([0 length(left_assrmi_noVR_meany)])
ylim([-0.25 0.25])
xlabel('sample');
ylabel('Value');
title('assrLI turning left');
grid on;
hold on
x_values = 1:length(left_assrmi_noVR_meany);
y_values = left_assrmi_noVR_mean;
plot(x_values, y_values, 'LineWidth', 2, 'Color', 'red');
highlight_patch = patch([lower_edge, upper_edge, upper_edge, lower_edge], [-1, -1, 1, 1], 'yellow');
set(highlight_patch, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off

subplot(2, 2, 2)
x_values = 1:length(right_assrmi_noVR_meany);
y_values = right_assrmi_noVR_meany;
plot(x_values, y_values,  'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
xlim([0 length(right_assrmi_noVR_meany)])
ylim([-0.25 0.25])
xlabel('sample');
ylabel('Value');
title('assrLI turning right');
grid on;
hold on
x_values = 1:length(right_assrmi_noVR_meany);
y_values = right_assrmi_noVR_mean;
plot(x_values, y_values, 'LineWidth', 2, 'Color', 'red');
highlight_patch = patch([lower_edge, upper_edge, upper_edge, lower_edge], [-1, -1, 1, 1], 'yellow');
set(highlight_patch, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off

% assrMI left zscore
subplot(2, 2, 3)
x_values = 1:length(left_assrmi_noVR_zscore);
y_values = left_assrmi_noVR_zscore;
plot(x_values, y_values,  'LineWidth', 2, 'Color', 'red');
xlim([0 length(left_assrmi_noVR_zscore)])
ylim([-2.5 2.5])
xlabel('sample');
ylabel('Value');
highlight_patch = patch([lower_edge, upper_edge, upper_edge, lower_edge], [-10 -10 10 10], 'yellow');
set(highlight_patch, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title('z-scored assrLI turning left');
grid on;

% assrMI right zscore
subplot(2, 2, 4)
x_values = 1:length(right_assrmi_noVR_zscore);
y_values = right_assrmi_noVR_zscore;
plot(x_values, y_values,  'LineWidth', 2, 'Color', 'red');
xlim([0 length(right_assrmi_noVR_zscore)])
ylim([-2.5 2.5])
xlabel('sample');
ylabel('Value');
highlight_patch = patch([lower_edge, upper_edge, upper_edge, lower_edge], [-10 -10 10 10], 'yellow');
set(highlight_patch, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title('z-scored assrLI turning right');
grid on;

tightfig;

% and save
output_file_path = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\new';
if ~exist(output_file_path, 'dir')      % Check if the folder exists, and create it if not
    mkdir(output_file_path);
end

savename = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\new\grand_mean_assrLI_new.png';
saveas(hFig, savename);

%% special new extra hot for thesis

lower_edge_sec = (lower_edge + length(assrPOW39_left_double)*0.2)*2;
upper_edge_sec = (upper_edge + length(assrPOW39_left_double)*0.2)*2;

mean_of_left = mean(left_assrmi_noVR_mean);
mean_of_right = mean(right_assrmi_noVR_mean);
left_mean_line = ones(1, length(left_assrmi_noVR_mean));
right_mean_line = ones(1, length(right_assrmi_noVR_mean));
left_mean_line = left_mean_line * mean_of_left;
right_mean_line = right_mean_line * mean_of_right;


hFig = figure; % Create a figure
subplot(1, 1, 1)

% Define the new x-values range
new_x_values = linspace(3045.2, 12180.8, length(left_assrmi_noVR_meany));

plot(new_x_values, left_assrmi_noVR_mean, 'LineWidth', 1, 'Color', 'blue');
hold on;
plot(new_x_values, left_mean_line, '--', 'LineWidth', 1, 'Color', 'blue');
xlim([min(new_x_values), max(new_x_values)]);
ylim([-0.1 0.1]);

% Set custom x-axis tick labels
xticks([new_x_values(1), mean(new_x_values), new_x_values(end)]);
xticklabels({'Curve Entry', 'Curve Apex', 'Curve Exit'});

% xlabel('Time (ms)');
ylabel('ASSR LI value');
% title('ASSR Lateralization Patterns');
grid on;


% Plot for right_assrmi_noVR
plot(new_x_values, right_assrmi_noVR_mean, 'LineWidth', 1, 'Color', 'red');
plot(new_x_values, right_mean_line, '--', 'LineWidth', 1, 'Color', 'red');
highlight_patch = patch([lower_edge_sec, upper_edge_sec, upper_edge_sec, lower_edge_sec], [-1, -1, 1, 1], 'yellow');
set(highlight_patch, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
legend('Turning left', 'Turning left mean', 'Turning right', 'Turning right mean', 'Location', 'northwest');
hold off

tightfig;

savename = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\new\grand_mean_assrLI_both.png';
saveas(hFig, savename);

%% get data for seperate powers

left_assr39_noVR = [];
left_assr41_noVR = []; 
right_assr39_noVR = []; 
right_assr41_noVR = [];

% data structure: data{sub, block}{left/right, 39/41}
for i = 1:length(assrMI_data_VR(:,1))
    for j = 1:length(assrMI_data_VR(1,:))
        left_assr39_noVR = [left_assr39_noVR; assrPOW_data_noVR{i,j}{1,1}];
        left_assr41_noVR = [left_assr41_noVR; assrPOW_data_noVR{i,j}{1,2}];
        right_assr39_noVR = [right_assr39_noVR; assrPOW_data_noVR{i,j}{2,1}];
        right_assr41_noVR = [right_assr41_noVR; assrPOW_data_noVR{i,j}{2,2}];
    end
end

left_assr39_noVR        = cell2mat(left_assr39_noVR);
left_assr41_noVR        = cell2mat(left_assr41_noVR);
right_assr39_noVR       = cell2mat(right_assr39_noVR);
right_assr41_noVR       = cell2mat(right_assr41_noVR);


%% take middle 60% to only look at curves 

index_of_20percent = round(length(left_assr39_noVR()) * 0.20); % same indices can be used for all data, since length should be the same after timewarping
index_of_80percent = round(length(left_assr39_noVR) * 0.80); % same indices can be used for all data, since length should be the same after timewarping

%% Prepare plot for seperate powers

% zscore
assrPOW39_left_zscore = zscore(mean(left_assr39_noVR, 1));
assrPOW39_right_zscore = zscore(mean(right_assr39_noVR, 1));
assrPOW41_left_zscore = zscore(mean(left_assr41_noVR, 1));
assrPOW41_right_zscore = zscore(mean(right_assr41_noVR, 1));

assrPOW39_left_zscore = assrPOW39_left_zscore(:, index_of_20percent:index_of_80percent);
assrPOW39_right_zscore = assrPOW39_right_zscore(:, index_of_20percent:index_of_80percent);
assrPOW41_left_zscore = assrPOW41_left_zscore(:, index_of_20percent:index_of_80percent);
assrPOW41_right_zscore = assrPOW41_right_zscore(:, index_of_20percent:index_of_80percent);

% Find the minimum value for each variable
min_values_left39 = min(assrPOW39_left_zscore);
min_values_left41= min(assrPOW41_left_zscore);
min_values_right39 = min(assrPOW39_right_zscore);
min_values_right41 = min(assrPOW41_right_zscore);

% Subtract the minimum value from each variable
assrPOW39_left_zscore = assrPOW39_left_zscore - min(min_values_left39);
assrPOW41_left_zscore = assrPOW41_left_zscore - min(min_values_left41);
assrPOW39_right_zscore = assrPOW39_right_zscore - min(min_values_right39);
assrPOW41_right_zscore = assrPOW41_right_zscore - min(min_values_right41);

%% stats for ASSR powers

% create struct for fieldtrip
datastr                     = [];
datastr.individual          = zeros(length(left_assr39_noVR(:,1)), 1, length(left_assr39_noVR));
datastr.time                = EEG_39_warped_left.times(:,index_of_20percent:index_of_80percent);
datastr.label               = {'N1'};
left39                      = datastr;
left39.individual(:,1,:)    = left_assr39_noVR;
left41                      = datastr;
left41.individual(:,1,:)    = left_assr41_noVR;

% run test
lengths                     = length(left_assr39_noVR(:,1));
cfg                         = [];
cfg.parameter               = 'individual';
cfg.latency                 = 'all'; %[-2 2];
cfg.method                  = 'montecarlo';
cfg.statistic               = 'depsamplesT';
cfg.alpha                   = 0.05;
cfg.correcttail             = 'prob';
cfg.numrandomization        = 1000;
cfg.correctm                = 'cluster';
cfg.design(1:2*lengths, 1)  = [ones(lengths, 1); 2*ones(lengths, 1)];
cfg.design(1:2*lengths, 2)  = [ones(lengths, 1); 1*ones(lengths, 1)];
cfg.ivar                    = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                    = 2; % the 2nd row in cfg.design contains the subject number
stat_left                   = ft_timelockstatistics(cfg, left39, left41);

%% ASSR powers separately 

hFig = figure; % Create a figure
screenSize = get(0, 'ScreenSize'); % Get the screen size
set(hFig, 'Position', screenSize);

% ASSR Powers left
subplot(2, 2, 1)
x_values_left = 1:length(mean(left_assr39_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_left = left_assr39_noVR(:, index_of_20percent:index_of_80percent);
plot(new_x_values, y_values_left, 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
xlim([min(new_x_values), max(new_x_values)]);
ylim([-0.2 0.2])
hold on;
x_values_left = 1:length(mean(left_assr39_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_left = mean(left_assr39_noVR(:, index_of_20percent:index_of_80percent), 1);
plot(new_x_values, y_values_left, 'LineWidth', 2, 'Color', 'red');
hold on;
x_values_left = 1:length(mean(left_assr39_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_left = -left_assr41_noVR(:, index_of_20percent:index_of_80percent);
plot(new_x_values, y_values_left, 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
hold on;
x_values_left = 1:length(mean(left_assr41_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_left = -mean(left_assr41_noVR(:, index_of_20percent:index_of_80percent), 1);
plot(new_x_values, y_values_left, 'LineWidth', 2, 'Color', 'red');
grid on;
title('ASSR powers turning left');
hold off;

% ASSR Powers right
subplot(2, 2, 2)
x_values_right = 1:length(mean(right_assr39_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_right = right_assr39_noVR(:, index_of_20percent:index_of_80percent);
plot(new_x_values, y_values_right, 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
xlim([min(new_x_values), max(new_x_values)]);
ylim([-0.2 0.2])
hold on;
x_values_right = 1:length(mean(right_assr39_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_right = mean(right_assr39_noVR(:, index_of_20percent:index_of_80percent), 1);
plot(new_x_values, y_values_right, 'LineWidth', 2, 'Color', 'red');
hold on;
x_values_right = 1:length(mean(right_assr39_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_right = -right_assr41_noVR(:, index_of_20percent:index_of_80percent);
plot(new_x_values, y_values_right, 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
hold on;
x_values_right = 1:length(mean(right_assr41_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_right = -mean(right_assr41_noVR(:, index_of_20percent:index_of_80percent), 1);
plot(new_x_values, y_values_right, 'LineWidth', 2, 'Color', 'red');
grid on;
title('ASSR powers turning right');
hold off;

% ASSR Powers left zscore
subplot(2, 2, 3)
x_values_left = 1:length(mean(assrPOW39_left_zscore, 1));
y_values_left = assrPOW39_left_zscore;
plot(new_x_values, y_values_left, 'LineWidth', 1, 'Color', 'red');
xlim([min(new_x_values), max(new_x_values)]);
hold on;
x_values_left = 1:length(mean(assrPOW39_left_zscore, 1));
y_values_left = -assrPOW41_left_zscore;
plot(new_x_values, y_values_left, 'LineWidth', 1, 'Color', 'red');
hold on;
grid on;
title('z-scored ASSR powers turning left');
hold off;

% ASSR Powers right zscore
subplot(2, 2, 4)
x_values_right = 1:length(mean(assrPOW39_right_zscore, 1));
y_values_right = assrPOW39_right_zscore;
plot(new_x_values, y_values_right, 'LineWidth', 1, 'Color', 'red');
xlim([min(new_x_values), max(new_x_values)]);
hold on;
x_values_right = 1:length(mean(assrPOW39_right_zscore, 1));
y_values_right = -assrPOW41_right_zscore;
plot(new_x_values, y_values_right, 'LineWidth', 1, 'Color', 'red');
hold on;
grid on;
title('z-scored ASSR powers turning right');
hold off;

tightfig;

% and save
output_file_path = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\new';
if ~exist(output_file_path, 'dir')      % Check if the folder exists, and create it if not
    mkdir(output_file_path);
end

savename = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\new\grand_mean_assrPOW_new.png';
saveas(hFig, savename);


%%


hFig = figure; % Create a figure
screenSize = get(0, 'ScreenSize'); % Get the screen size
set(hFig, 'Position', screenSize);

% ASSR Powers left
subplot(2, 2, 1)
x_values_left = 1:length(mean(left_assr39_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_left = mean(left_assr39_noVR(:, index_of_20percent:index_of_80percent), 1);
plot(x_values_left, y_values_left, 'LineWidth', 2, 'Color', 'red');
xlim([0 length(mean(left_assr39_noVR(:, index_of_20percent:index_of_80percent), 1))])
ylim([-0.06 0.06])
hold on;
x_values_left = 1:length(mean(left_assr41_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_left = -mean(left_assr41_noVR(:, index_of_20percent:index_of_80percent), 1);
plot(x_values_left, y_values_left, 'LineWidth', 2, 'Color', 'red');
hold on;
grid on;
title('ASSR powers turning left');
hold off;

% ASSR Powers right
subplot(2, 2, 2)
x_values_right = 1:length(mean(right_assr39_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_right = right_assr39_noVR(:, index_of_20percent:index_of_80percent);
plot(x_values_right, y_values_right, 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
xlim([0 length(mean(right_assr39_noVR(:, index_of_20percent:index_of_80percent), 1))])
ylim([-0.2 0.2])
hold on;
x_values_right = 1:length(mean(right_assr39_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_right = mean(right_assr39_noVR(:, index_of_20percent:index_of_80percent), 1);
plot(x_values_right, y_values_right, 'LineWidth', 2, 'Color', 'red');
hold on;
x_values_right = 1:length(mean(right_assr39_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_right = -right_assr41_noVR(:, index_of_20percent:index_of_80percent);
plot(x_values_right, y_values_right, 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
hold on;
x_values_right = 1:length(mean(right_assr41_noVR(:, index_of_20percent:index_of_80percent), 1));
y_values_right = -mean(right_assr41_noVR(:, index_of_20percent:index_of_80percent), 1);
plot(x_values_right, y_values_right, 'LineWidth', 2, 'Color', 'red');
grid on;
title('ASSR powers turning right');
hold off;

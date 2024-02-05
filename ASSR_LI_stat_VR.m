%% first run PIPELINE_assrLI_all

%% get data

left_assrmi_VR = [];
right_assrmi_VR = []; 
left_assrmi_noVR = []; 
right_assrmi_noVR = [];

% data structure: data{sub, block}{left/right}
for i = 1:length(assrMI_data_VR(:,1))
    for j = 1:length(assrMI_data_VR(1,:))
        left_assrmi_VR = [left_assrmi_VR; assrMI_data_VR{i,j}{1,1}];
        right_assrmi_VR = [right_assrmi_VR; assrMI_data_VR{i,j}{1,2}];
        left_assrmi_noVR = [left_assrmi_noVR; assrMI_data_noVR{i,j}{1,1}];
        right_assrmi_noVR = [right_assrmi_noVR; assrMI_data_noVR{i,j}{1,2}];
    end
end

%% cut

left_assrmi_VR = left_assrmi_VR(:,index_of_20percent:index_of_80percent);
right_assrmi_VR = right_assrmi_VR(:,index_of_20percent:index_of_80percent);

left_assrmi_VR_mean = mean(left_assrmi_VR, 1);
right_assrmi_VR_mean = mean(right_assrmi_VR, 1);
left_assrmi_VR_zscore = zscore(left_assrmi_VR_mean);
right_assrmi_VR_zscore = zscore(right_assrmi_VR_mean);

%% statistics: t-test for VR

left_assrmi_VR_mean2              = mean(left_assrmi_VR, 2);
right_assrmi_VR_mean2             = mean(right_assrmi_VR, 2);
right_assrmi_VR_mean2_fitting     = right_assrmi_VR_mean2(1:length(left_assrmi_VR_mean2),:); % match size
left_mean_mean                      = mean(left_assrmi_VR_mean2);
right_mean_mean                     = mean(right_assrmi_VR_mean2);

left                                = left_assrmi_VR_mean2;
right                               = right_assrmi_VR_mean2_fitting;

[h, p, ci, stats] = ttest(left, right, 'tail', 'right', 'alpha', 0.05);
ttest_VR.h = h; 
ttest_VR.p = p;
ttest_VR.ci = ci;
ttest_VR.stats = stats;

if p > 0.05
    fprintf('The difference between left and right assrMI values did not reach statistical significance (p = %f)!!\n', ttest_VR.p);
else
    fprintf('The difference between left and right assrMI values did reach statistical significance (p = %f)!!\n', ttest_VR.p);
end

mean_left_ttest = mean(left_assrmi_VR_mean2);
std_left_ttest = std(left_assrmi_VR_mean2);
mean_right_ttest = mean(right_assrmi_VR_mean2_fitting);
std_right_ttest = std(right_assrmi_VR_mean2_fitting);


%% statistics: permutation VR

% create struct for fieldtrip
datastr                 = [];
datastr.individual      = zeros(length(left_assrmi_VR(1:150,1)), 1, length(left_assrmi_VR));
datastr.time            = EEG_39_warped_left.times(:,index_of_20percent:index_of_80percent);
datastr.label           = {'N1'};
% datastr.dimord        = 'subj_chan_time';
left                    = datastr;
left.individual(:,1,:)  = left_assrmi_VR(1:150,:);
right                   = datastr;
right_assrmi_VR_fitting = right_assrmi_VR(1:length(left_assrmi_VR(:,1)),:);
right.individual(:,1,:) = right_assrmi_VR_fitting(1:150,:);

% run test
lengths                 = length(left_assrmi_VR(1:150,1));
cfg                     = [];
cfg.parameter           = 'individual';
cfg.latency             = 'all'; %[-2 2];
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT';
cfg.alpha               = 0.05;
cfg.correcttail         = 'prob';
cfg.numrandomization    = 1000;
cfg.correctm            = 'cluster';
cfg.design(1:2*lengths, 1)   = [ones(lengths, 1); 2*ones(lengths, 1)];
cfg.design(1:2*lengths, 2)   = repelem(1:6, 1, 50)';
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
stat_VR               = ft_timelockstatistics(cfg, left, right);



%% plot assrMI grand mean

hFig = figure; % Create a figure
screenSize = get(0, 'ScreenSize'); % Get the screen size
set(hFig, 'Position', screenSize);

% assrMI left
subplot(3, 2, 1)
x_values = 1:length(left_assrmi_VR);
y_values = left_assrmi_VR;
plot(x_values, y_values,  'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
xlim([0 length(left_assrmi_VR)])
ylim([-1 1])
xlabel('sample');
ylabel('Value');
title('assrLI turning left');
grid on;
hold on
x_values = 1:length(left_assrmi_VR);
y_values = left_assrmi_VR_mean;
plot(x_values, y_values, 'LineWidth', 2, 'Color', 'red');
highlight_patch = patch([lower_edge, upper_edge, upper_edge, lower_edge], [-1, -1, 1, 1], 'yellow');
set(highlight_patch, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off

subplot(3, 2, 2)
x_values = 1:length(right_assrmi_VR);
y_values = right_assrmi_VR;
plot(x_values, y_values,  'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
xlim([0 length(right_assrmi_VR)])
ylim([-1 1])
xlabel('sample');
ylabel('Value');
title('assrLI turning right');
grid on;
hold on
x_values = 1:length(right_assrmi_VR);
y_values = right_assrmi_VR_mean;
plot(x_values, y_values, 'LineWidth', 2, 'Color', 'red');
highlight_patch = patch([lower_edge, upper_edge, upper_edge, lower_edge], [-1, -1, 1, 1], 'yellow');
set(highlight_patch, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off

% assrMI left zscore
subplot(3, 2, 3)
x_values = 1:length(left_assrmi_VR_zscore);
y_values = left_assrmi_VR_zscore;
plot(x_values, y_values,  'LineWidth', 2, 'Color', 'red');
xlim([0 length(left_assrmi_VR_zscore)])
ylim([-2.5 2.5])
xlabel('sample');
ylabel('Value');
highlight_patch = patch([lower_edge, upper_edge, upper_edge, lower_edge], [-10 -10 10 10], 'yellow');
set(highlight_patch, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title('z-scored assrLI turning left');
grid on;

% assrMI right zscore
subplot(3, 2, 4)
x_values = 1:length(right_assrmi_VR_zscore);
y_values = right_assrmi_VR_zscore;
plot(x_values, y_values,  'LineWidth', 2, 'Color', 'red');
xlim([0 length(right_assrmi_VR_zscore)])
ylim([-2.5 2.5])
xlabel('sample');
ylabel('Value');
highlight_patch = patch([lower_edge, upper_edge, upper_edge, lower_edge], [-10 -10 10 10], 'yellow');
set(highlight_patch, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title('z-scored assrLI turning right');
grid on;

% check normal distribution
subplot(3, 2, 5)
histogram(left_assrmi_VR, 'Normalization', 'probability', 'EdgeColor', 'w');
title('Histogram turning left');

subplot(3, 2, 6)
histogram(right_assrmi_VR, 'Normalization', 'probability', 'EdgeColor', 'w');
title('Histogram turning right');

tightfig;

% and save
output_file_path = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\new';
if ~exist(output_file_path, 'dir')      % Check if the folder exists, and create it if not
    mkdir(output_file_path);
end

savename = 'C:\Users\ericw\Desktop\Master_Thesis\Data\CURRENT_DATA\4_FIGURES\new\grand_mean_assrLI_new.png';
saveas(hFig, savename);
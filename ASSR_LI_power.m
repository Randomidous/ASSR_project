%% get data for seperate powers

left_assr39_noVR = [];
left_assr41_noVR = []; 
right_assr39_noVR = []; 
right_assr41_noVR = [];

% data structure: data{sub, block}{left/right, 39/41}
for i = 1:length(assrMI_data_VR(:,1))
    for j = 1:length(assrMI_data_VR(1,:))
        left_assr39_noVR = [left_assr39_noVR; mean(cell2mat(assrPOW_data_noVR{i,j}{1,1}),1)];
        left_assr41_noVR = [left_assr41_noVR; mean(cell2mat(assrPOW_data_noVR{i,j}{1,2}),1)];
        right_assr39_noVR = [right_assr39_noVR; mean(cell2mat(assrPOW_data_noVR{i,j}{2,1}),1)];
        right_assr41_noVR = [right_assr41_noVR; mean(cell2mat(assrPOW_data_noVR{i,j}{2,2}),1)];
    end
end

%% take middle 60% to only look at curves 

index_of_20percent = round(length(left_assr39_noVR()) * 0.20); % same indices can be used for all data, since length should be the same after timewarping
index_of_80percent = round(length(left_assr39_noVR) * 0.80); % same indices can be used for all data, since length should be the same after timewarping

left_assr39_noVR = left_assr39_noVR(:,index_of_20percent:index_of_80percent);
left_assr41_noVR = left_assr41_noVR(:,index_of_20percent:index_of_80percent);
right_assr39_noVR = right_assr39_noVR(:,index_of_20percent:index_of_80percent);
right_assr41_noVR = right_assr41_noVR(:,index_of_20percent:index_of_80percent);

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
right39                     = datastr;
right39.individual(:,1,:)   = right_assr39_noVR;
right41                     = datastr;
right41.individual(:,1,:)   = right_assr41_noVR;

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
cfg.design(1:2*lengths, 2)  = [repelem(1:6, 1, 1)'; repelem(1:6, 1, 1)'];
cfg.design                  = cfg.design';
cfg.ivar                    = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                    = 2; % the 2nd row in cfg.design contains the subject number

stat_left                   = ft_timelockstatistics(cfg, left39, left41);
stat_right                  = ft_timelockstatistics(cfg, right39, right41);

%% Prepare plot for seperate powers

% zscore
assrPOW39_left_zscore = zscore(mean(left_assr39_noVR, 1));
assrPOW39_right_zscore = zscore(mean(right_assr39_noVR, 1));
assrPOW41_left_zscore = zscore(mean(left_assr41_noVR, 1));
assrPOW41_right_zscore = zscore(mean(right_assr41_noVR, 1));

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


%% ASSR powers separately 

hFig = figure; % Create a figure
screenSize = get(0, 'ScreenSize'); % Get the screen size
set(hFig, 'Position', screenSize);

% Define the new x-values range
new_x_values = linspace(3045.2, 12180.8, length(left_assr39_noVR));

% Define the custom x-axis tick positions and labels
custom_xticks = [min(new_x_values), mean(new_x_values), max(new_x_values)];
custom_xticklabels = {'Curve Entry', 'Curve Apex', 'Curve Exit'};

% ASSR Powers left
% subplot(2, 2, 1)
figure; 
x_values_left = 1:length(mean(left_assr39_noVR, 1));
y_values_left = left_assr39_noVR;
plot(new_x_values, y_values_left,'--', 'LineWidth', 1, 'Color', [0.4, 0.4, 1]);
hold on;
x_values_left = 1:length(mean(left_assr39_noVR , 1));
y_values_left = mean(left_assr39_noVR , 1);
plot(new_x_values, y_values_left, 'LineWidth', 2, 'Color', 'blue');
hold on;
x_values_left = 1:length(mean(left_assr39_noVR , 1));
y_values_left = -left_assr41_noVR ;
plot(new_x_values, y_values_left, '--', 'LineWidth', 1, 'Color', [1, 0.4, 0.4]);
hold on;
x_values_left = 1:length(mean(left_assr41_noVR , 1));
y_values_left = -mean(left_assr41_noVR , 1);
plot(new_x_values, y_values_left, 'LineWidth', 2, 'Color', 'red');
grid on;
title('ASSR powers turning left');
% Set custom x-axis tick labels
xticks(custom_xticks);
xticklabels(custom_xticklabels);
xlim([min(new_x_values), max(new_x_values)]);
ylim([-0.15 0.15]);
ylabel('µV²/Hz');


% ASSR Powers right
figure;
% subplot(2, 2, 2)
x_values_right = 1:length(mean(right_assr39_noVR, 1));
y_values_right = right_assr39_noVR;
plot(new_x_values, y_values_right, '--', 'LineWidth', 1, 'Color', [0.4, 0.4, 1]);
hold on;
x_values_right = 1:length(mean(right_assr39_noVR , 1));
y_values_right = mean(right_assr39_noVR , 1);
plot(new_x_values, y_values_right,  'LineWidth', 2, 'Color', 'blue');
hold on;
x_values_right = 1:length(mean(right_assr39_noVR , 1));
y_values_right = -right_assr41_noVR ;
plot(new_x_values, y_values_right,  '--', 'LineWidth', 1, 'Color', [1, 0.4, 0.4]);
hold on;
x_values_right = 1:length(mean(right_assr41_noVR , 1));
y_values_right = -mean(right_assr41_noVR , 1);
plot(new_x_values, y_values_right, 'LineWidth', 2, 'Color', 'red');
grid on;
title('ASSR powers turning right');
% Set custom x-axis tick labels
xticks(custom_xticks);
xticklabels(custom_xticklabels);
xlim([min(new_x_values), max(new_x_values)]);
ylim([-0.05 0.05]);
ylabel('µV²/Hz');

% ASSR Powers left zscore
figure;
% subplot(2, 2, 3)
x_values_left = 1:length(mean(assrPOW39_left_zscore, 1));
y_values_left = assrPOW39_left_zscore;
plot(new_x_values, y_values_left, 'LineWidth', 1, 'Color', 'blue');
% xlim([min(new_x_values), max(new_x_values)]);
hold on;
x_values_left = 1:length(mean(assrPOW39_left_zscore, 1));
y_values_left = -assrPOW41_left_zscore;
plot(new_x_values, y_values_left, 'LineWidth', 1, 'Color', 'red');
hold on;
grid on;
title('z-scored ASSR powers turning left');
xlim([min(new_x_values), max(new_x_values)]);
xticks(custom_xticks);
xticklabels(custom_xticklabels);
hold off;


% ASSR Powers right zscore
figure;
% subplot(2, 2, 4)
x_values_right = 1:length(mean(assrPOW39_right_zscore, 1));
y_values_right = assrPOW39_right_zscore;
plot(new_x_values, y_values_right, 'LineWidth', 1, 'Color', 'blue');
% xlim([min(new_x_values), max(new_x_values)]);
hold on;
x_values_right = 1:length(mean(assrPOW39_right_zscore, 1));
y_values_right = -assrPOW41_right_zscore;
plot(new_x_values, y_values_right, 'LineWidth', 1, 'Color', 'red');
hold on;
grid on;
title('z-scored ASSR powers turning right');
xlim([min(new_x_values), max(new_x_values)]);
xticks(custom_xticks);
xticklabels(custom_xticklabels);
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
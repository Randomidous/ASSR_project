%% export fooof data for ANOVA in R

%% Here we leave matlab and go to a better place

variablesToKeep = {'fooof_results_noVR', 'fooof_results_VR'};
allVariables = who;
variablesToDelete = setdiff(allVariables, variablesToKeep);
for i = 1:length(variablesToDelete)
    eval(['clear ' variablesToDelete{i}]);
end

goodSubs                = {1, 2, 8, 9, 10, 11};
VR                      = {1, 2, 9};
noVR                    = {8, 10, 11};
condition               = {'_A_VR','_B_VR','_A','_B'};

%% Create table for ANOVA

% Create a data table

n = length(goodSubs);           % the amount of participants
nGroups = 2;                    % Vr and noVR
nDirection = 2;                 % left and right (surprise!)
nCond = length(condition);     % won't comment

participants = repelem(1:n, nDirection * nCond)'; % Two measurements per participant
groups = repelem({'Vr', 'no_VR'}, 1, n * nDirection)';
blocks = repmat(repelem(1:nCond, 1, nDirection), 1, n)';
directions = repmat({'L', 'R'}', n * nCond, 1);

power_cell = cell(1, length(goodSubs));

for isub = 1:length(goodSubs)
    for iblock = 1:length(condition)
        power_cell{1,isub} = cell(1, 4);
        power_cell{1,isub}{1,1} = power_total{isub,iblock}{1,1};
        power_cell{1,isub}{1,2} = power_total{isub,iblock}{1,2};
        power_cell{1,isub}{1,3} = power_total{isub,iblock}{1,1}; % NOT CORRECT
        power_cell{1,isub}{1,4} = power_total{isub,iblock}{1,2}; % NOT CORRECT
    end
end
power = cat(2, power_cell{:,:})';
power = cell2mat(power);

dataTable = table(participants, groups, blocks, directions, power);

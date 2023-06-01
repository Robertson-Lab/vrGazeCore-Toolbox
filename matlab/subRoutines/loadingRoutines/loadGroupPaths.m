groupHeatPlotsDir = [paths.projectHeatPlotsDir params.cohortName '/'];
groupHeatMatDir = [paths.projectHeatMatDir params.cohortName '/'];
%timecourse
% groupTimecourseMatDir = [paths.projectTimecourseMatDir params.cohortName '/']; %timecourse mat files
% groupTimecoursePlotsDir = [paths.projectTimecoursePlotsDir params.cohortName '/']; %timecourse plots


% list of paths
groupPathsList = {
    'groupHeatPlotsDir'; 
    'groupHeatMatDir';
    % 'groupTimecourseMatDir';
    % 'groupTimecoursePlotsDir';
    }; 

%% Add group paths to struct for easy passing into functions
for pathI = 1:length(groupPathsList) % loop through all paths
    paths.(char( groupPathsList(pathI) ) )=eval(char( groupPathsList(pathI) ));
end

%% Check if necessary directories exist; if not, make them.
    for pathI = 1:length(groupPathsList) % loop through all paths
        if ~exist( eval(char( groupPathsList(pathI) )) ) % if path does not exist
            mkdir( eval(char( groupPathsList(pathI) )) ) % create it!
        end
    end
clearvars -except paths params
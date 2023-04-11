function [paths] = loadSubjectPaths(inputSubject,paths,params)

%fixation directory
subjectFixMatDir = [paths.projectFixMatDir inputSubject '/'];
subjectFixPlotsDir = [paths.projectFixPlotsDir inputSubject '/'];
%saccade directory
%subjectSaccadeDir = [paths.projectSaccadeDataDir inputSubject '/'];
%heatmaps
subjectHeatMatDir = [paths.projectHeatMatDir inputSubject '/']; %heatmap mat files
subjectHeatPlotsDir = [paths.projectHeatPlotsDir inputSubject '/']; %heatmap plots
%timecourse
% subjectTimecourseMatDir = [paths.projectTimecourseMatDir inputSubject '/']; %timecourse mat files
% subjectTimecoursePlotsDir = [paths.projectTimecoursePlotsDir inputSubject '/']; %timecourse plots

subPathsList = {
    'subjectFixMatDir';
    'subjectFixPlotsDir'
    %'subjectSaccadeDir';
    'subjectHeatMatDir';
    'subjectHeatPlotsDir';
    % 'subjectTimecourseMatDir';
    % 'subjectTimecoursePlotsDir';
};

%% Add group paths to struct for easy passing into functions
for pathI = 1:length(subPathsList) % loop through all paths
    paths.(char( subPathsList(pathI) ) )=eval(char( subPathsList(pathI) ));
end

%% Check if necessary directories exist; if not, make them.
    for pathI = 1:length(subPathsList) % loop through all paths
        if ~exist( eval(char( subPathsList(pathI) )) ) % if path does not exist
            mkdir( eval(char( subPathsList(pathI) )) ) % create it!
        end
    end

clearvars -except paths params
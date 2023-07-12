%% Project Parameters Load File %%
% DESCRIPTION:
% This file sets paths to be passed to other scripts and adds dirs to your path.

% Set Top Level Paths
projectDir = '~/Documents/matGazeCore/';                       % point this to to your main project folder (one level above folders with stimuli, raw data, etc.)
gazeCoreDir = [projectDir 'matlab/'];                       % point this to the vrGazeCore script directory 

% Set stimuli and data folder to use                    
projectRawDataDir = [projectDir 'rawData/'];                    % Raw Data
projectStimDir = [ projectDir 'stimuli/combined/' ];            % Project stimuli to process

%% Based on the top-level folders, populate path name variables:
% Set analysis results folders
projectDataDir = [ projectDir 'eyeTrackResults/' ];                               
    projectFixDataDir = [ projectDataDir 'fixations/' ];         % Fixation .mat files
        projectFixMatDir = [ projectFixDataDir 'mat/' ]; 
        projectFixPlotsDir = [ projectFixDataDir 'plots/' ];
    %projectSaccadeDataDir = [ projectDataDir 'saccades/' ];             % Saccade info
    projectHeatDir = [ projectDataDir 'heatMaps/' ];                      % Heat maps
        projectHeatMatDir = [ projectHeatDir 'mat/' ]; 
        projectHeatPlotsDir = [ projectHeatDir 'plots/' ]; 
            projectTimecoursePlotsDir = [ projectDataDir 'timecourseHeat/' ]; % Timecourse gifs of heat maps

% Set meta and log folders
projectAnalLogsDir = [ projectDir 'eyeTrackLogs/' ];
    projectMetaDataDir = [ projectAnalLogsDir 'meta/' ];               % Meta Data
    projectLogsDir = [ projectAnalLogsDir 'logs/' ];                   % Logs

%% Print check about paths that need to be changed
fprintf('Check that the following paths are correct:\n'); % Display necessary paths
    fprintf('projectDir = %s\n',projectDir);
    fprintf('gazeCoreDir = %s\n',gazeCoreDir);
    fprintf('projectRawDataDir = %s\n',projectRawDataDir);
    fprintf('projectStimDir = %s\n',projectStimDir);
% If correct, input 1
checkPath = input('Are the paths correct? \n 1 = Yes \n 2 = No\n Enter:');
if checkPath == 2
    % break script, prompt to go back
    fprintf('batchProcessPars will stop running. Go to setPaths to change the incorrect paths.\n')
    error('Paths are not correct.')
end

% list of paths
pathsList = {
    'projectDir'; 
    'gazeCoreDir'; 
    'projectRawDataDir';
    'projectStimDir'; 
    'projectDataDir'; 
    'projectFixDataDir'; 
    'projectFixMatDir';
    'projectFixPlotsDir';
 %   'projectSaccadeDataDir';
    'projectHeatDir';
    'projectHeatMatDir';
    'projectHeatPlotsDir';
    'projectTimecoursePlotsDir';
    'projectAnalLogsDir';
    'projectMetaDataDir'; 
    'projectLogsDir';
    }; 

%% Store paths in a struct for easy passing into functions
paths = struct;
for pathI = 1:length(pathsList) % loop through all paths
    paths.(char( pathsList(pathI) ) )=eval(char( pathsList(pathI) ));
end


%% Check if necessary directories exist; if not, make them.
    for pathI = 1:length(pathsList) % loop through all paths
        if ~exist( eval(char( pathsList(pathI) )) ) % if path does not exist
            mkdir( eval(char( pathsList(pathI) )) ) % create it!
        end
    end


%% add important things to  MATLAB path
addpath(genpath(paths.gazeCoreDir)) % Scripts
addpath(genpath(paths.projectDataDir)) % Project data
addpath(genpath(paths.projectRawDataDir)) % Raw Data
addpath(genpath(paths.projectStimDir))  % Stimulus files
addpath(genpath(paths.projectAnalLogsDir)) % Project data


%% Clear vars
clearvars -except paths %Clear everything not in a struct
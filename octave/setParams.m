%% Project Parameters File %%
% DESCRIPTION:
% This file contains all of the parameters for fixation and heat analysis.
% Modify these values here. They will be passed on to the other scripts via
% the params and paths structs.

%% Release Version
scriptVersion = '0.1.5'; % don't change unless updating version

%% Need to Adjust Parameters

    % Processing Options
    runFindFix = 1; % if 1, runs find fixations for all individuals
    runHeatmappingIndivid = 1; % if 1, run individual heatmapping for each subject
    runHeatmappingGroup = 1; % if 1, run group heatmapping for the cohort
    runTimecourseGifIndivid = 0; % if 1, make individual heatmapping timecourse gifs for each subject
    runTimecourseGifGroup = 0; % if 1, make individual heatmapping timecourse gifs for each subject

    % Headset Type
    headsetType = 0; %Default 0; DK2=0 Vive=1 ViveEye = 2, Oculus Go = 3

    % Select Subjects
    cohortName = 'octTest-1'; % pick a cohort name. If doing group processing, this will be used to label your files
    listSubjectNames = 0; % if 1 list

    if listSubjectNames == 1 % manually list subjects here
        subjectNames = cellstr({
            'subject1';
        });
    else % OR
        %%%% Run All Subjects in a specified raw Data Dir
        d = dir(sprintf(paths.projectRawDataDir));%list directory where live
        subjectFiles = d(~ismember({d.name}, {'.', '..','.DS_Store'}));%exclude non-files
        subjectNames = [];
        for subjectI = 1:length(subjectFiles)
            subjectName = subjectFiles(subjectI).name(1:end-4);
            subjectNames = [subjectNames ; cellstr(subjectName)];%record the current scene to the list of scenes
        end
    end

    % Heatmapping Options
    heatmapTimesteps = [1]; % how many time segments to divide scene into

    % Scene Parameter
    sceneLength = 16; % scene duration in seconds

%% Additional Scene Parameters
    % Scene Sample Filter
    minSamples = 100; %Default 100; skip scenes w/ less than 100 samples
    % Scene Time Filter
    minTimeInScene = 5; %Default 5; skip scenes shorter than 5 seconds.
    % scannedFilter
    scannedFilter = 1; % if set to 1, removes scenes where they were not explored
    scannedThresh = 66; % default 66, set to a percent threshold to exclude scenes not explored by that much based on head direction

    % Manually Exclude Scenes %
    % Each subject gets a string of scenes to excluded. Subjects are seperated by a comma
    excludeScenes = { '','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','' };
        %%%%%%%%%  EXAMPLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % excludeScenes = { 'scene01 scene02 scene05' , ...                                   %%
        %     '', '', 'scene05'};                                                             %%
        % this means par1 exclude 3 scenes. pars2+3 exclude none. par4 exclude only scene05   %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ignoreList = {'fixate*'}; % anything in scene directory that should be ignored for all subjects (e.g., subfolders)

    % Load scenes based on the sceneDir (make sceneList)
    d = dir(sprintf(paths.projectStimDir));%list directory where live
    sceneFiles = d(~ismember({d.name}, {'.', '..','allScenes.mat', '.DS_Store', ignoreList{:} }));%exclude non-files
    sceneList = [];
    for sceneI = 1:length(sceneFiles)
        sceneName = sceneFiles(sceneI).name(1:end-4);
        sceneList = [sceneList ; cellstr(sceneName)];%record the current scene to the list of scenes
    end

    junkSceneList = {'endTrial'}; % need to check if this can be merged with ignoreList

%% Eyetracking options
        gazeType = 0; %Default 0 (2D tracking), 3D = 1
    % Fixations Options
        excludeFirstNSec = 0; % default = 0; exclude n seconds of start of trial
        minMad = 50; % default = 50; Windows with a MAD (mean absolute deviation) less than x deg/s were classified as potential fixations
        maxMad = 100; % default = 100;  Windows with a MAD more than x deg/s were classified as likely saccades
        excludeFixDursLessThan = 0.1; % exclude any fixation with a duration less than 0.1s
    % Eye Selection
        useEye = 3; % if stereo eye tracking, 0=eye0, 1=eye1, 2=choose best eye at each point, 3 = average eyes and set confidence to lowest confidence of 2 eyes
    % Locked Head to center?
        headLocked = 0; % default 0. If this is 1, sets all head direction to straight ahead. Do this if your data was head locked
    % Confidence Filter
        minConfThresh = 0.25; %Default 0.25; Data with lower than this pupil confidence value will be discarded by confidenceFilter; default for old pupil should be around 0.6 but with new pupil version, lower like 0.25
        maxConfPercent = 75; %Default 75; If more than 75% of the data is discarded due to confidence, skip the scene
    % Eccentricity Filter
        eccFiltX = 100; %set the max eccentricty in the x dimension that will be analyzed (gaze points falling farther than +/- eccFiltX/2 will be discarded)
        eccFiltY = 100; %gaze points falling farther than +/- eccFiltY/2 on the y axis will be discarded
                            % in the old code this was 100x100 and corresponded to roughly dva. in the
                            % new scene viewing it is in proportion to headset "fov" value so for dk2 it is now
                            % proportional to ~ 145x160 ? confusing
    % Smoothing
        useSmoothing = 0; % default 0. if 1, use Matt's smoothing. Not exactly fully implemented at the moment but could be easily. [DP: figure out where smoothing is]
        useInterpolation = 0; % default 0, if 1 use linear interpolation (Wass 2012) to recover data loss samples < 100ms
        durationForInterpolation = 0.10; %using 100ms
    % Fixation Calculation
        % based on Mean Absolute Deviation in degrees/sec
        fixType = 1; % type of fixation to calculate, 1 = gaze fixation, 2 = head fixation
            if fixType == 1 % for eye fixations
                fixSpatialDist = 2; % if fixations within this # of deg, group
                fixTempDist = 0.15; % if fixations are within fixTempDist of each other group
            elseif fixType == 2
                fixSpatialDist = 2; % if fixations within this # of deg, group
                fixTempDist = 0.10; % if fixations are within fixTempDist of each other, group
            end
        fixValidation = 0; % 1 = validate fixations with saccades on each side
        valWindow = 0; %number of MAD samples to average both preceding and following a given fixation for validation
        excludeLastNSec = 0; % default = 0; exclude n seconds of end of trial
        excludeFirstNFix = 1; % default = 1; exclude n fixations at start of trial
        excludeLastNFix = 0; % default = 0; exclude n fixations at end of trial
        valWindow = 2; %number of points before a fixation to check
    % Drift Correction
        driftCorrection = 0; % if set to 1, then do drift correction based on last pretrial
    % Pre-Trial Fixations
        concatSanity = 1; %Set Equator scenes to the correct name for plotting! we do actually want to do this don't remove
        avgPreTrial = 1; % If set to 1, calculate pretrial fixation distance and threshold
        avgPreTrialThresh = 5; % DVA, tolerance for determining if a scene should be drift corrected(belowthresh) or skipped (above thresh), distance from true point (drift)
        excludeByPreTrial = 1; % exclude the next scene if the previous pre-trial was bad
        avgPreTrialX = -1.26; % true x coordinate of stimulus -1.26
        avgPreTrialY = -0.54; % true y coordinate of stim -0.54

%% Plotting & Heatmap Options
    % General Options
        plotVisibility = 'off'; %if 'off', figures are not visible when making them
    % Fixation Plotting
        plotHeadFixationInfo = 0; % Default 0; plot head yaw path (w/ fix points), velocity, acceleration, MAD on image
        plotCombinedGazeFlag = 0; % Default 1; plot combined head, gaze, fix points on image
        plotHeadRawFlag = 0; % Default 0; plot the raw HMD direction center point on image
        plotGazeRawFlag = 0; % Default 0; plot the raw gaze points on image
        plotFixFlag = 1; % Default 1; plot the fixations on image
        plotMADFlag = 0; % Default 1; plot viewport position & MAD for each trial over time
    % Heatmap Options
        trimFactor = 0; % if zero, don't trim, else trim top and bottom by trim factor
        padAmount = 100; % Default = 100; set the amount of padding around the border ( to account for smoothing at edges )
        boundFiltering = 1; % Default 1; if it is 0, don't do upper and lower bound filtering as in group method. If 1, do it.
    % GIF options
        delayTime = 0.1; % delay time between images in gif
        deleteTimecourseJPEG = 0; % if 1, delete timecourse .jpeg images after making gif
    % Equirectangulage image parameters
        imDimX = 2000; %horizontal dimension (in pixels) of equirectangular images
        imDimY = 1000; %vertical dimension (in pixels) of equirectangular images

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fixed Parameters
%% Probably don't modify; change at your own risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Headset-Dependent Settings

if headsetType==2 %Vive Eye settings
    avgPreTrialX = -1.26; % true x coordinate of stimulus -1.26
    avgPreTrialY = -0.54; % true y coordinate of stim -0.54
    if gazeType == 1
        driftCorrection = 0;
        avgPreTrial = 1;
    end
end

if headsetType == 3 %if we're using the oculus go, override any other settings here
    fixType = 2; % make sure we're using head fixations
    useSmoothing = 0; %make sure no smoothing
    useInterpolation = 0; % make sure no interpoolation
    avgPreTrial = 0;

    % don't actually know these numbers - need to measure
    fovX = 100;%horizontal fov of headset in degrees
    fovY = 100;%vertical fov of headset in degrees
end

%%  Set headset effective "fov" based on code version and headset type JSM CHECK THIS
if headsetType == 0 % DK2
    %measured fov numbers
    fovX = 90;%horizontal "fov" of headset in degrees
    fovY = 100;%vertical "fov" of headset in degrees

%         % jeff expansion fov numbers
%         fovX = 145;%horizontal "fov" of headset in degrees
%         fovY = 160;%vertical "fov" of headset in degrees
elseif headsetType == 1 % VIVE
    %note, these are probably not correct. If you use the Vive for some reason, double check this...
    fovX = 145;%horizontal "fov" of headset in degrees 4/30//19 was 145
    fovY = 160;%vertical "fov" of headset in degrees4/30//19 was 160

elseif headsetType == 2 % VIVE EYE
    fovX = 75;%horizontal "fov" of headset in degrees
    fovY = 90;%vertical "fov" of headset in degrees
end

if fovX > fovY
    maxFOV = fovX;
else
    maxFOV = fovY;
end


%% Parameter Check
% Print paths to be checked
fprintf('Check that the following parameters are correct:\n'); % Display necessary paths
    fprintf('runFindFix = %d\n',runFindFix);
    fprintf('runHeatmappingIndivid = %d\n',runHeatmappingIndivid);
    fprintf('runHeatmappingGroup = %d\n',runHeatmappingGroup);
    fprintf('runTimecourseGifIndivid = %d\n',runTimecourseGifIndivid);
    fprintf('heatsetType = %d\n',headsetType);
    fprintf('cohortName = %s\n',cohortName);
    fprintf('subjectNames = %s\n',subjectNames{:})
    fprintf('gazeType = %d\n',gazeType);
% If correct, input 1
checkParams = input('Are the parameters correct? \n 1 = Yes \n 2 = No\n Enter:');
if checkParams == 2
    % error, prompt to edit file
    fprintf('batchProcessPars will stop running. Go to setParams to change the incorrect parameters.\n')
    error('Parameters are not correct.')
end


%% List of parameters in this file:
% if you add any new parameters to this file, also add them to this list so
% that they are added to the params struct
paramsList = {
    'scriptVersion';
    'cohortName';
    'subjectNames';
    'excludeScenes';
    'runFindFix';
    'runHeatmappingIndivid';
    'runHeatmappingGroup';
    'headsetType';
    'gazeType';
    'useEye';
    'sceneLength';
    'minSamples';
    'minTimeInScene';
    'minConfThresh';
    'maxConfPercent';
    'eccFiltX';
    'eccFiltY';
    'concatSanity';
    'avgPreTrialThresh';
    'useSmoothing';
    'useInterpolation';
    'durationForInterpolation';
    'fixType';
    'fixValidation';
    'valWindow';
    'minMad';
    'maxMad';
    'excludeFirstNSec';
    'excludeLastNSec';
    'excludeFirstNFix';
    'excludeLastNFix';
    'excludeFixDursLessThan';
    'avgPreTrial';
    'avgPreTrialX';
    'avgPreTrialY';
    'junkSceneList';
    'fixSpatialDist';
    'fixTempDist';
    'valWindow';
    'plotHeadFixationInfo';
    'plotCombinedGazeFlag';
    'plotHeadRawFlag';
    'plotGazeRawFlag';
    'plotMADFlag';
    'padAmount';
    'boundFiltering';
    'sceneFiles';
    'sceneList';
    'plotFixFlag';
    'fovX';
    'fovY';
    'maxFOV';
    'imDimX';
    'imDimY';
    'excludeByPreTrial';
    'headLocked';
    'driftCorrection';
    'scannedFilter';
    'scannedThresh';
    'trimFactor';
    'heatmapTimesteps';
    'plotVisibility';
    'runTimecourseGifIndivid';
    'runTimecourseGifGroup';
    'delayTime';
    'deleteTimecourseJPEG';
    };

% Store params in a struct for easy passing into functions
params = struct;
for paramI = 1:length(paramsList) % loop through all paths
    params.(char( paramsList(paramI) ) )=eval(char( paramsList(paramI) ));
end
params = orderfields(params); %alphabetize structure

clearvars -except paths params %Clear everything not in a struct

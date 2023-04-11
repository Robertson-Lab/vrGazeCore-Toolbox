% Function to go from calculated fixations and durations to individual and or group 
% heatmap .mat files and plots [as in Henderson & Hayes (2018)]
% INSTRUCTIONS:
% This script can be called from the batchProcessPars.m script
% findFixations needs to be run before this step.
function fix2Heat(inputSubjects,paths,params)


%% Create Log File, log things
fileID = fopen([ paths.projectLogsDir, 'fix2heat', datestr(now,'mmmmdd-yyyy_HH-MM-SS'),'.txt'],'w');%create log file
fprintf(fileID, '\nvrGazeCore Version: %s', params.scriptVersion); %print script version to log
fprintf(fileID, '\nStart Time: %s\nParams:\n', datetime); %print start datetime to log
paramsList = fieldnames(params);
for iParam = 1:length(paramsList) % log all of the params
    if isstruct( params.(char(paramsList(iParam))) )
        continue % skip the param if it is a struct
    end
    fprintf(fileID, '%s: %s\n', string(paramsList(iParam)), string( params.(char(paramsList(iParam))) ) ); %print script version to log
end
clear paramsList iParam
%% Detect if there was a list or a single subject as input, then do group or individ processing
if size(inputSubjects,1)==1
    groupHeat = 0;
else
    groupHeat = 1;
end
fprintf(fileID, '\ngroupHeat = %s', string(groupHeat)); %print script version to log
%% Process data
% for each scene, find the upper and lower bounds to scale each particpants data
% add all subject gaze, create one duration weighted group heatmap
for sceneIndex=1:length(params.sceneList) % cycle through scenes
    % sceneDurs = []; % create variable to hold within scene duration data
    % sceneX = []; % create variable to hold within scene x fixation position
    % sceneY = [];% create variable to hold within scene y fixation position
    currentScene = char( params.sceneList(sceneIndex) ); %get current scene name
    currentSceneStripped = currentScene(1:end-3); %strip the scene order number from end of name DEEPA: is this necessary
    fprintf(fileID,'\n\n*~*~*~*~* Processing segment: %03d\n', sceneIndex);%log current scene# to console
    fprintf('\n\n*~*~*~*~* Processing segment: %03d\n', sceneIndex);%log current scene# to console
    fprintf(fileID,'Scene: %s\n', currentScene);%log current scene name to console
    fprintf('Scene: %s\n', currentScene);%log current scene name to console

    %% Aggregate fixations across subjects for each timestep
    % go through each subject & divide data into specified number of timesteps

    % need to add an extra to get bounding range 
    timestepBounds = linspace(0, params.sceneLength, params.heatmapTimesteps+1);

    % timestep array will be time (rows) x items (1=fix x, 2=fix y, 3=durations)
    timestepArray = cell(params.heatmapTimesteps, 1);

    for parIndex = 1:size(inputSubjects,1) % cyle through participants
        
        subjectName = char (inputSubjects(parIndex,:) ); % get subject name

        % load data for current subject --> we do this here to avoid
        % loading multiple times
        [subjectFixationStartTimes, subjectFixationEquiX, subjectFixationEquiY, subjectFixationDurations] = getParticipantFixations(subjectName, currentScene, fileID, paths, params);


        % add subject data to time matrix for current scene
        for step = 1:params.heatmapTimesteps
            % find the fixations that are greater than the current start
            % and less than the next start time
            timestepIdxs = subjectFixationStartTimes >= timestepBounds(step) & subjectFixationStartTimes < timestepBounds(step+1);

            % get the current fixations and durations for the timestep
            subjectFixationEquiX(timestepIdxs);
            subjectFixationEquiY(timestepIdxs);
            stepFixationDurations = subjectFixationDurations(timestepIdxs);

            % we bound filter and scale the durations within the current timestep
            stepFixationDurations = scaleFixationDurations(stepFixationDurations, params);

            % create a matrix for the current timestep
            stepFixationInformation = [
                subjectFixationEquiX', subjectFixationEquiY', stepFixationDurations'
            ];

            % append to the timestep array
            timestepArray{step} = cat(1, timestepArray{step}, stepFixationInformation);
        end
    end
    
    %% Create heatmaps
    % for each timestep create a heatmap

    fprintf(fileID,'\n\n----Heatmapping----\n');%log current scene# to console
    fprintf('\n\n----Heatmapping----\n');%log current scene# to console

    timestepHeatmaps = cell(params.heatmapTimesteps, 1);

    for step = 1:params.heatmapTimesteps
        if isempty(timestepArray{1}) % skip if this subject didn't see this scene and if no subjects saw this scene and doing group maps DEEPA check if should be step == 1
            % add array of zeros here since there isnt a map
            timestepHeatmaps{step} = cat(1, timestepHeatmaps{step}, nan(2000, 1000));
            fprintf(fileID,'\nStep ', step, '/', params.heatmapTimesteps, ' Empty. Skipping.');%log current scene# to console
            continue
        end

        % if the step exists, create a heatmap for the timestep
        % timestep array will be time (rows) x items (1=fix x, 2=fix y, 3=durations)
        sceneHeat = createSceneHeatmap(timestepArray{step}(:,1), timestepArray{step}(:,2), timestepArray{step}(:,3), params.trimFactor);
        
        % append to timestep heatmap array
        timestepHeatmaps{step} = cat(1, timestepHeatmaps{step}, sceneHeat);
    end

    % Save timestep heatmap array
    if groupHeat == 1
        saveName = [paths.groupHeatMatDir '/' 'groupheat_' params.cohortName '_timestep' num2str(params.heatmapTimesteps) '-' num2str(step) '_' currentScene '.mat'];
    else
        saveName = [paths.subjectHeatMatDir subjectName '_timestep' num2str(params.heatmapTimesteps) '-' num2str(step) '_' currentScene '.mat'];
    end
    save(saveName,'timestepHeatmaps');

    %% Plot heatmaps

    imFullName = dir([paths.projectStimDir currentScene '*.jp*']);
    imFullName = imFullName.name;

    % skip if the image is not there
    if isempty(imFullName)
        continue
    end

    if isnan(timestepHeatmaps{step})
        continue
    end

    sceneImage = imread([paths.projectStimDir imFullName]);
    sceneImage = imresize(sceneImage,[1000,2000]); % resize scene image

    for step = 1:params.heatmapTimesteps
        % plotting takes in the matrix of scene heat
        % for each time step plots the image 
        fig = figure('visible',params.plotVisibility), imshow(sceneImage); %show with fixed aspect ratio,etc
        hold on

        stepHeatmap = pcolor(timestepHeatmaps{step});

        shading flat;
        colormap(parula(500));
        alpha = (~isnan(timestepHeatmaps{step}))*0.2;
        hold on
        set(stepHeatmap,'facealpha',0.7);
        ax=gca;
        ax.Visible = 'off';%no borders
        f=getframe;%get just image
        
        if groupHeat == 1
            saveName = [ paths.groupHeatPlotsDir '/' 'groupheat_' params.cohortName '_timestep' num2str(params.heatmapTimesteps) '-' num2str(step) '_' currentScene '.jpeg'];
        else
            saveName = [ paths.subjectHeatPlotsDir '/' subjectName '_timestep' num2str(params.heatmapTimesteps) '-' num2str(step) '_' currentScene '.jpeg'];
        end

        imwrite(f.cdata,saveName);%write image to jpeg
        close all

    end

    clearvars -except fileID params paths sceneIndex parIndex groupHeat inputSubjects
end
fprintf(fileID, '\n\nEnd Time: %s\n', datetime);
end
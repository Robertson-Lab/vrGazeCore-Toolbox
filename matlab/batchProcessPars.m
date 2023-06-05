%% Batch Processing Participants %%
% DESCRIPTION:
% Script to BATCH process participants through fixation calculation and fixation density mapping
% Loads participant names from setParams and runs through all of them

%% Set Up Directories
setPaths
setParams

%% Process subject data
for parIndex = 1:size(params.subjectNames,1)
    % set up subject
    subjectName = params.subjectNames(parIndex,:); % get current subject name
    excludedScenesPar = string(params.excludeScenes(parIndex)); % get excluded scenes list for current subject
    paths = loadSubjectPaths(char(subjectName),paths,params);

    %% Run findFixations (Individual)
    if params.runFindFix == 1
        findFixations(char(subjectName),paths,params) % run findFixation for current subject
        fprintf('Finished finding fixations for subject %s\n', subjectName)
    end
    
    %% Run Heatmapping (Individual)
    if params.runHeatmappingIndivid == 1 %[DP: and runHeatmappingGroup == 0]
        fix2Heat(subjectName,paths,params)
        fprintf('Finished generating heatmaps for subject %s\n', subjectName)
    end

    % Run gif making
    if params.runTimecourseGifIndivid == 1
        plotTimecourseGif(char(subjectName), paths,params)
        fprintf('Finished making GIFs for subject %s\n', subjectName)
    end
    
    fprintf('Finished processing subject %s\n', subjectName)


    clearvars -except parIndex paths params
end

%% Process group data
if (params.runHeatmappingGroup == 1 || params.runTimecourseGifGroup)
    loadGroupPaths

    % Run Heatmapping (Group)
    if params.runHeatmappingGroup == 1
        fix2Heat(params.subjectNames,paths,params)
        fprintf('Finished generating heatmaps for cohort %s\n', params.cohortName)
    end

    %Run gif making (Group)
    if params.runTimecourseGifGroup == 1
        plotTimecourseGif(params.cohortName,paths,params)
        fprintf('Finished making GIFs for cohort %s\n', params.cohortName)
    end

    fprintf('Finished processing cohort %s\n', params.cohortName)
end



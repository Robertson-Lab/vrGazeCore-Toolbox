%% DP's Notes:
%  PURPOSE: to define fixations from VR eye-tracking data
%  INPUTS: 1) participant name; 2) excluded scene strings for participant; 3) raw data file from VR eye-tracker output
%  OUTPUTS: TBD 

function findFixations(subjectName,paths,params)
%script to define fixations following Peterson (2016), Wass (2012)
% INPUTS: 1) participant name; 2) excluded scenes string for that participant
% loads raw data file from pupil labs output as struct ds: txt file with 8-10 columns:
% [sceneName, date, 
% timestamp (relative to comp), timesetamp (relative to trialstart, 
% eye0x (normalized coords, 0-1), eye0x (normalized coords, 0-1), 
% eye1x (normalized coords, 0-1), eye1y (normalized coords, 0-1),
% eye0 confidence, eye1 confidence]

%% Create Log File, Log Everything
uniqueID = [subjectName, datestr(now,'mmmmdd-yyyy_HH-MM-SS')];
fileID = fopen([ paths.projectLogsDir, uniqueID, 'findfix.txt'],'w');%create log file
fprintf(fileID, '\nvrGazeCore Version: %s', params.scriptVersion); %print script version to log
fprintf(fileID, '\nStart Time: %s\nParams:\n', datetime); %print start datetime to log

%% Create META File, log METADATA things
metaData = struct; % instantiate structure for meta data
paramsList = fieldnames(params);
for iParam = 1:length(paramsList) % log all of the params
    if isstruct( params.(char(paramsList(iParam))) )
        continue % skip the param if it is a struct though
    end
    fprintf(fileID, '%s: %s\n', string(paramsList(iParam)), string( params.(char(paramsList(iParam))) ) ); %print script version to log
end

%% Load raw data, organize into a matrix, separate data into scenes

[parData, sceneNames, sceneChangeList] = processRawData(subjectName, params, paths);

%% start variables to hold things
skipNextScene = 0;
timeList = [];
fpsList = [];
sceneProcessedList = [""];
sceneCount = 0;
lowConfidenceCount = 0;
excludedPretrialCount = 0;
eyeTrackerFailCount = 0;
excludedScannedCount = 0;
minSampleCount = 0;
%notScannedCount = 0;
percentScannedList = [];
meanRawConfidenceList =  [];

junkSceneList = {'cloudland', 'recalibrate', 'endTrial'};
% set drift correction amounts to 0 so there is at least a value for the
% first scene
preshiftX = 0;
preshiftY = 0;

bom = '﻿';

preFixCount = 0; %counter for pretrial fixation scenes
%% loop through scenes
for s=1:length(sceneChangeList)-1 %loop through scenes
    
    currentScene = sceneNames{sceneChangeList(s)+1};
    currentSceneText = [currentScene sprintf('%03d',s) ];
    sceneProcessedList = [sceneProcessedList; (currentScene);];
    fprintf(fileID,'\n\n*~*~*~*~* Processing segment: %03d\n', s);%log current scene# to console
    fprintf('\n\n*~*~*~*~* Processing segment: %03d\n', s);%log current scene# to console
    fprintf(fileID,'Scene: %s\n', currentScene);%log current scene name to console
    fprintf('Scene: %s\n', currentScene);%log current scene name to console
    
    if strcmp(currentScene, '_sanityTarget360') || strcmp(currentScene, 'fixate') % are we in a pre-trial scene
        preTrialScene = 1; % turn on preTriaScene flag then
        preFixCount = preFixCount + 1;
    else
        preTrialScene = 0;% otherwise turn off the preTrialScene flag
        sceneCount = sceneCount+1;
    end        
    
    % trim parData to only considered scenes
    % the new unity gaze logging is different than the old version. Account for this.
    if params.unityProjectVersion == 0 %old unity. 16p,activepass,pilots
        parDataScene = parData(sceneChangeList(s):sceneChangeList(s+1),:);
    elseif params.unityProjectVersion == 1 %new template tlb, emax,snapi,newer
        parDataScene = parData((sceneChangeList(s)+1):sceneChangeList(s+1),:);
    end
    
    
    timeInScene = parDataScene(end,1)-parDataScene(1,1); %find scene length (seconds)
    metaData.(['scene' currentSceneText]).timeInScene = timeInScene;
    samplesInScene = length(parDataScene);
    metaData.(['scene' currentSceneText]).samplesInScene = samplesInScene;
    
    %calculate scene fps - col 1 is time
    sceneFPS = mean(1./diff(parDataScene(:,1)));
    fpsList(end+1) = sceneFPS; %tack on the current fps to the list
    timeList(end+1) = timeInScene; %tack on the current time to the list
    metaData.(['scene' currentSceneText]).samplesInScene = samplesInScene;

%   skip manually excluded scenes (specified by excludedScenesPar in coreParams) [DP check this]
%   junk scenes (specified line 46 in findFixations --> related to old data collection)
%   and scenes less than minTimeInScene
    if any(strcmp(currentScene,junkSceneList) | timeInScene < params.minTimeInScene) % TLB need to figure out strcmp(currentScene, params.excludedScenesPar) || 
       continue
    end
    
    if skipNextScene == 1  && params.excludeByPreTrial == 1 && preTrialScene==0 %if there was a bad pre-trial scene last time & we are excludingbypretrial & currently not a pretrialscene
        fprintf(fileID,'Scene Excluded by bad Pre-Trial! Skipping.\n');
        fprintf('Scene Excluded by bad Pre-Trial! Skipping.\n');
        metaData.(['scene' currentSceneText]).excludedByPreTrial = 1;
        metaData.(['scene' currentSceneText]).included = 0;
        excludedPretrialCount = excludedPretrialCount+1;
        continue %skip
    end
    
    fprintf(fileID,'Scene Duration: %.2f seconds\n', timeInScene);%log time in scene
    fprintf('Scene Duration: %.2f seconds\n', timeInScene);%log time in scene
    fprintf(fileID,'Scene FPS: %.2f frames per second\n', sceneFPS);%log fps
    fprintf('Scene FPS: %.2f frames per second\n', sceneFPS);%log fps
    [rawSceneData] = parseRaw2struct(parDataScene,params.fovX,params.fovY,preTrialScene,params,preshiftX,preshiftY); %tlb 4-23-19 changing parData to parDataScene to only output scenedata 
    if params.scannedFilter == 1
        if preTrialScene == 0
            yawScan = parDataScene(:,3);
            yawScan = yawScan + 180;
            percentScanned = 100 * ( max(yawScan)-min(yawScan) ) / 360;
            fprintf(fileID,'Scanned %.2f percent of scene\n', percentScanned);%log fps
            fprintf('Scanned %.2f percent of scene\n', percentScanned);%log fps
            percentScannedList(end+1) = percentScanned;
            if percentScanned < params.scannedThresh
                excludedScannedCount = excludedScannedCount+1;
                continue
            end
            %notScannedCount = notScannedCount+1;
        end
    end
    
    %%%%%%%% Discard Beginning and End Seconds  %%%%%%%%%%
    %77 frames per second.. (MIT Desktop +  DK2 + old unity code...)
    % exclude   start +  fps * seconds          to       end       - fps *    seconds
    parDataScene = parDataScene(1+( round(sceneFPS)*params.excludeFirstNSec) : length(parDataScene) - ( round(sceneFPS)* params.excludeLastNSec), :);


    if params.headsetType ~= 3 %don't filter data for oculus go
        
        %%%%%%%% TLB TESTING 2D Pre-Trial based Drift Correction 
        % if we are not in a pre-trial and we are doing drift correction, then do it
        if preTrialScene == 0 && params.driftCorrection == 1
            parDataScene(:,5) = -preshiftX + parDataScene(:,5);
            parDataScene(:,6) = -preshiftY + parDataScene(:,6);
        end

        %%%%%%%% CONFIDENCE FILTER %%%%%%%%%%
        if preTrialScene == 0 % if it is not a pre trial scene, log the mean raw confidence
            meanRawConfidenceList(end+1) = mean( parDataScene(:,7) );
        end

        %tlb 8-29-19 - slightly changing filtering method to keep all data but label invalid data
        %filter out values below confidence threshold

        [parDataScene,confidenceFlag,confidenceRemoved,confidencePercent] = confidenceFilter(parDataScene,params,fileID);


        metaData.(['scene' currentSceneText]).confidenceRemoved = confidenceRemoved;
        metaData.(['scene' currentSceneText]).confidencePercent = confidencePercent;

        if confidenceFlag == 1 % if removing more than 75% then skip that scene
            fprintf(fileID,'MORE THAN %d PERCENT REMOVED! Skipping.\n',params.maxConfPercent);
            fprintf('MORE THAN %d PERCENT REMOVED! Skipping.\n',params.maxConfPercent);
            metaData.(['scene' currentSceneText]).excludedConfidence = 1; 
            metaData.(['scene' currentSceneText]).included = 0;
            lowConfidenceCount = lowConfidenceCount + 1;
            if preTrialScene == 1  && params.excludeByPreTrial == 1 % if we are throwing out the pretrial scene for low confidence, then skpi the next real scene
                skipNextScene = 1;
            end
            continue %skip
        else
            metaData.(['scene' currentSceneText]).excludedConfidence = 0;
        end

        %%%%%%%% ECCENTRICITY FILTER %%%%%%%%%%
        [parDataScene,ecc,percentRemoved] = eccFilter(parDataScene,params.eccFiltX,params.eccFiltY,params.fovX,params.fovY,fileID);
        eccentricityRemoved = ecc;
        percentRemoved;

        %tlb 8-29-2019 - added invalid index column to raw scene data
        filteredData = parDataScene;
        filteredData(filteredData(:,end) == 1, 1:end-1) = nan; %
    %     filteredData = parDataScene(parDataScene(:,end) ~= 1,1:end-1);  % then remove rows labelled invalid by filters
    
        %tlb just changing out variable names
        if length(filteredData) < params.minSamples
            fprintf(fileID,'NOT ENOUGH SAMPLES! Skipping\n');
            fprintf('NOT ENOUGH SAMPLES! Skipping\n');
            metaData.(['scene' currentSceneText]).excludedEccentricity = 1; 
            metaData.(['scene' currentSceneText]).included = 0;

            if preTrialScene == 1 % if we are throwing out the pretrial scene for no samples, then skip the next real scene
                skipNextScene = 1;
            else % if not a pretrial, increase excluded count list
                minSampleCount = minSampleCount+1;
            end
            continue %skip
        else
            metaData.(['scene' currentSceneText]).excludedConfidence = 0;
        end

        if length(unique(filteredData(:,6))) < 10 % if there are less than 10 unique gaze positions then something broke!
            fprintf(fileID,'EYE TRACKER FAILED! Skipping\n');
            fprintf('EYE TRACKER FAILED! Skipping\n');
            metaData.(['scene' currentSceneText]).excludedEyeTrackerFailed = 1;
            metaData.(['scene' currentSceneText]).included = 0;
            if preTrialScene == 1  && params.excludeByPreTrial == 1 % if we are throwing out the pretrial scene, then skip the next real scene
                skipNextScene = 1;
            else
                eyeTrackerFailCount = eyeTrackerFailCount+1;
            end
            continue % NeXt!
        else
            metaData.(['scene' currentSceneText]).excludedEyeTrackerFailed = 0;
        end

    else
        
        filteredData = parDataScene;
        
    end
    

    %%%%%%%% PARSE DATA INTO VARIABLES %%%%%%%%
    [time,pitch,yaw,roll_rad,gazeX_deg,gazeY_deg,gaze_ecc,confidence] = parseDS(filteredData,params);

    
    %%%%%%%% Head Locked Option %%%%%%%%
    if params.headLocked == 1 || strcmp(currentScene, 'fixate')
        pitch(:) = 90; % if we are head locked, set all head direction to straight ahead
        yaw(:) = 180;
        roll_rad(:) = 0;
    end
    
    %%%%%%%% RECTIFY/ROTATE %%%%%%%% 
    % This is where eye coordinates are convolved with head coordinates to
    % produce gaze on a sphere as pitch and yaw
    if params.headsetType ~= 3
        [gaze_yaw_sphere,gaze_pitch_sphere] = rectifyGaze(gazeX_deg, roll_rad, gazeY_deg, pitch, yaw); % converts the gaze points from viewport to sphere space

        % jsm fix this
        rawSceneData.gaze_yaw_sphere = gaze_yaw_sphere - 180;
        rawSceneData.gaze_pitch_sphere = gaze_pitch_sphere - 90;


    % % % % % % % % % % % % % %        % ARM EDIT HERE TO HACKILY change the manipulations that happen to
    % % % % % % % % % % % % % %    % gazepitchsphere and gazerawsphere
    % % % % % % % % % % % % % %     if params.fixType == 2 
    %         rawSceneData.rawHeadYaw = rawSceneData.rawHeadYaw + 180;
    %         rawSceneData.rawHeadPitch = rawSceneData.rawHeadPitch + 90;

    % % % % % % % % % % % % % %         rawSceneData.gaze_pitch_sphere = rawSceneData.rawHeadPitch;
    % % % % % % % % % % % % % %     end 



        %%%%%%%% WRAP POINTS %%%%%%%% TLB - can we just use wrapTo360?
        [gaze_yaw_sphere,gaze_pitch_sphere] = wrapPointsEquirect( gaze_yaw_sphere, gaze_pitch_sphere, 360, 180, 0); % don't round
    end
    
    %%%%%%%%%% BILATERAL SMOOTHING ALA PETERSON 2016 %%%%%%%%
    %using a bilateral filtering algorithm which eliminates jitter in
    %fixation while preserving saccades excludes samples where only one eye
    %is available data is converted into a single x and y position averaged
    %across both eyes source: http://people.csail.mit.edu/jiawen/#code need
    %to figure out the parameters for this function and if there need to be
    %any tweaks
    % output = bilateralFilter( data, edge, sigmaSpatial, sigmaRange, ...
    %     samplingSpatial, samplingRange )
    
    if params.useSmoothing == 1
        numSamples = 3;
        
        % we bilateral filter w/ nans (which are removed samples by the
        % confidence/eccentricity filter) to prevent temporally distant
        % samples from being smoothed together (which occurs when data is
        % removed instead of labelled nan).

        [J, abs_vel, abs_acc, n] = bilatFilt(gaze_yaw_sphere,gaze_pitch_sphere,time,numSamples,numSamples,numSamples,0);
        
        %tlb updates 8-29-2019
        %remove nans from beginning and end of arrays that are a result of
        %the filter being bilateral (requires numSamples before and after a
        %given smoothing point).
        
        gaze_yaw_sphere = J.x(numSamples+1:end-numSamples)'; %assign smoothed values to vars used by calculate fixations
        gaze_pitch_sphere = J.y(numSamples+1:end-numSamples)';
        
        %adjust time and other fixation data to match length of gaze coordinate array
        time = time(numSamples+1:end-numSamples);
        gazeX_deg = gazeX_deg(numSamples+1:end-numSamples);
        gazeY_deg = gazeY_deg(numSamples+1:end-numSamples);
        gaze_ecc = gaze_ecc(numSamples+1:end-numSamples);
        confidence = confidence(numSamples+1:end-numSamples);
        
    end

    %%%%%%%%%  INTERPOLATION  %%%%%%%%%%%%
    % Written by TLB 2019
    
    if params.useInterpolation == 1
        
        validIdxs = [];
        
        %preserve the invalid samples (nans) as a reference array for
        %interpolating
        preInterpYaw = gaze_yaw_sphere;
        preInterpPitch = gaze_pitch_sphere;
        
        % remove invalid samples (nans) from all fixation related data in preparation for
        % interpolating between invalid samples
        gaze_yaw_sphere = gaze_yaw_sphere(~isnan(gaze_yaw_sphere));
        gaze_pitch_sphere = gaze_pitch_sphere(~isnan(gaze_pitch_sphere));
        
        time = time(~isnan(time));
        gazeX_deg = gazeX_deg(~isnan(gazeX_deg));
        gazeY_deg = gazeY_deg(~isnan(gazeY_deg));
        gaze_ecc = gaze_ecc(~isnan(gaze_ecc));
        confidence = confidence(~isnan(confidence));
        
        % these idxs refer where samples were removed from the raw data
        removedIdxs = find(parDataScene(:,end) == 1);

        
        %find number of idxs removed previous to a point
        try 
            if ~isempty(removedIdxs)
                if params.useSmoothing == 0 %if no bilateralFilter was used, we have no need to shift idxs
                    numSamples = 0;
                end
                
                %shift the rawTime (used to find the missing data b/c var
                % "time" is filtered) by numSamples removed from beginning by bilateralFilter
                alignRawTime = rawSceneData.rawTime(numSamples+1:end-numSamples);

                idxDiff = [1 ; diff(removedIdxs)];
                
                %find groups of missing data by looking for differences of removed indices >1
                %then shift indices by the number of samples removed from the beginning by bilateralFilter
                %to align the removed idxs with the filtered data
                startInterp = [removedIdxs(find(idxDiff > 1))] - numSamples; 
                startInterp = [removedIdxs(1) - numSamples ; startInterp(1:end-1)];
                
                %repeat for finding the end of missing data --> look for
                %the idx before a separation of missing data >1. shift for
                %bilateral filter
                endInterp = [removedIdxs(find(idxDiff > 1) - 1)] - numSamples; 

                if isempty(endInterp) %if there is only one removed sample, endInterp will be empty
                    endInterp = [removedIdxs(end) - numSamples];
                end

                lenInterp = []; %keep track of the length of each interpolation --> used to increment the idxs for interpolation

                % in the second column, add 0's to denote valid samples
                validIdxs = repmat(0, size(gaze_yaw_sphere,1), 1); % ones will be added where interpolations occured

                for idx = 1:size(startInterp,1) % go through each possible interpolation

                    if alignRawTime(endInterp(idx)) - alignRawTime(startInterp(idx)) < params.durationForInterpolation % if the amount of time is < 100ms 
                        try
                            % find point before and point after data was removed
                            alignToRawIdxs = [alignRawTime(startInterp(idx)-1:endInterp(idx)+1)];

                            startInterpIdx = find(ismember(time,alignToRawIdxs(1))); %find the point before the removed data in the processed data
                            endInterpIdx = find(ismember(time,alignToRawIdxs(end))); % point after the removed data (will be the consecutive point in the filtered data)

                            %because the data we are looking at is filtered, we only need the start index
                            interpX = [gaze_yaw_sphere(startInterpIdx:endInterpIdx)]';
                            interpY = [gaze_pitch_sphere(startInterpIdx:endInterpIdx)]';

                            %time of sample before removed data, time of sample after removed data
                            surroundTimes = [time(startInterpIdx);time(endInterpIdx)];

                            interpX(isnan(interpX)) = [];
                            interpY(isnan(interpY)) = [];

                            interpXY = interp1(surroundTimes, [interpX',interpY'], alignToRawIdxs);

                            %trim the values interpolated removing the reference points of valid data
                            interpXY = interpXY(2:end-1,:);

                            %then insert values into matrix of data
                            validIdxs = [ validIdxs(1:startInterp(idx) - 1) ; repmat(1, size(interpXY,1),1) ; validIdxs(startInterp(idx):end)]; %mark where the data was inserted with a 1
                            gaze_yaw_sphere = [ gaze_yaw_sphere(1:startInterp(idx) - 1) ; interpXY(:,1) ; gaze_yaw_sphere(startInterp(idx):end)];
                            gaze_pitch_sphere =  [ gaze_pitch_sphere(1:startInterp(idx) - 1) ; interpXY(:,2) ; gaze_pitch_sphere(startInterp(idx):end)];
                            time =  [ time(1:startInterp(idx) - 1) ; alignRawTime(startInterp(idx):endInterp(idx)) ; time(startInterp(idx):end)];

                            lenInterp(idx) = (endInterp(idx) + 1) - startInterp(idx); % use to increment interpolation idxs after adding values
                        catch
                            
                        end
                    end
                end   
            end
        catch
            
        end
    else
        removedIdxs = [];
        
        if params.headsetType ~= 3
            % remove invalid samples (nans) from all fixation related data 
            gaze_yaw_sphere = gaze_yaw_sphere(~isnan(gaze_yaw_sphere));
            gaze_pitch_sphere = gaze_pitch_sphere(~isnan(gaze_pitch_sphere));

            time = time(~isnan(time));
            gazeX_deg = gazeX_deg(~isnan(gazeX_deg));
            gazeY_deg = gazeY_deg(~isnan(gazeY_deg));
            gaze_ecc = gaze_ecc(~isnan(gaze_ecc));
            confidence = confidence(~isnan(confidence));
            
        end
    end
    
    %%%%%%%%% CALCULATE FIXATIONS FOR SPHERE DATA  %%%%%%%%%
    if params.fixType == 1 % gaze data 
        [mean_fix_yaw,mean_fix_pitch,fix_duration,length_fix,start_time,end_time,begin_fix,end_fix,mean_fix_gaze_ecc] = calculateFixations(removedIdxs, gazeX_deg,gazeY_deg,gaze_yaw_sphere,gaze_pitch_sphere,time,gaze_ecc,fileID,paths,subjectName,currentSceneText,params,rawSceneData);
    elseif params.fixType == 2 % head data 
         [mean_fix_yaw,mean_fix_pitch,fix_duration,length_fix,start_time,end_time,begin_fix,end_fix,mean_fix_gaze_ecc] = calculateFixations(removedIdxs, gazeX_deg,gazeY_deg,rawSceneData.rawHeadYaw,rawSceneData.rawHeadPitch,rawSceneData.rawTime,gaze_ecc,fileID,paths,subjectName,currentSceneText,params,rawSceneData);
    end
    
    %%% ARM MAKE SURE THAT Deg doesnt get used 
    %%%% figure out what is valid indexes 
    %%%% use raw time 
    
   %%% use pitch and yaw for more complicated way; pitch and yaw are
   %%% already filtered by confidence 
    
    
    %%%%%%%%% CALCULATE SACCADES ON VIEWPORT SPACE (VPS) DATA %%%%%%%%%
    % maybe put this inside calculateFixations so not to have unused vars in memory
%     if preTrialScene == 0 % only want to find saccades for gaze movement and non-sanity scenes
%         [sacc_durs,length_sacc,sacc_start_times, sacc_end_times,begin_sacc_idx,end_sacc_idx,mean_sacc_gaze_ecc] = calculateSaccades(gazeX_deg,gazeY_deg,time,gaze_ecc,fileID,paths,subjectName,currentSceneText,params,rawSceneData);
%     end
    
    %%Pretrial fixation error calculation
    if params.avgPreTrial == 1
        if preTrialScene == 1
            %if we want to calculate pretrial, we're not on the first
            %sanity scene (equator scene), and the current scene is a
            %sanity scene
% 
            nSec = 0.05;
            
            if params.gazeType == 1 %because Vive Eye coords are already in degrees, we want to use the raw deg coordinates
                preFixX = rawSceneData.rawViewportX_deg;
                preFixY = rawSceneData.rawViewportY_deg;  
            else
                preFixX = rawSceneData.rawViewportX_norm;
                preFixY = rawSceneData.rawViewportY_norm;      
            end

            preTime = rawSceneData.rawTime;
            
            endIdx = find(abs(preTime - (preTime(end) - nSec)) == min(abs(preTime - (preTime(end) - nSec)))); 
            
            
            for point = 1:size(preFixX,1)
                currEnd = find(abs(preTime - (preTime(point) + nSec)) == min(abs(preTime - (preTime(point) + nSec)))); 
                
                if currEnd >= endIdx
                    continue;
                else
                    currSampleX = mean(preFixX(point:currEnd));
                    currSampleY = mean(preFixY(point:currEnd));
                    alignTime(point) = preTime(point);
                    [fixShift(point), foo] = distance(currSampleX, currSampleY, 0, 0); 
                    
                    %adding to shift points on 2d screen space
                    fixShiftX(point) = currSampleX;
                    fixShiftY(point) = currSampleY;
                end
            end
            
            startIdx = find(abs(alignTime - (alignTime(end) - 5)) == min(abs(alignTime - (alignTime(end) - 5))));
            
            fixShift = mean(fixShift(startIdx:end));
            fixShiftX = mean(fixShiftX(startIdx:end));
            fixShiftY = mean(fixShiftY(startIdx:end));
            
%             fixShiftX = viewportNorm2dva(fixShiftX);
%             fixShiftY = viewportNorm2dva(fixShiftY);
%             fixShift = distance(fixShiftX, fixShiftY, params.avgPreTrialX, params.avgPreTrialY);
            
            fixShift = distance(fixShiftX, fixShiftY, 0, 0);
            
            if params.gazeType == 1 %distance is already in dva
                fixShift = fixShift;
            else
                fixShift = viewportNorm2dva(fixShift);
            end
            
%             % tlb - i think we should change the way to 
%             fixShift = rad2deg(2*atan2(fixShift,0.8));
            
            xShiftList(preFixCount,1) = fixShiftX;
            yShiftList(preFixCount,1) = fixShiftY;
            fixShift =  mean(fixShift);
            fixShiftList(preFixCount,1) = fixShift;
            if fixShift > params.avgPreTrialThresh && preTrialScene == 1 && params.excludeByPreTrial == 1 % if we are throwing out the pretrial scene, then skip the next real scene
                skipNextScene = 1;
                fprintf('PRE-TRIAL FIXATION: mean distance from center was over threshold: %4f', fixShiftList(preFixCount,1));
                fprintf('SKIP NEXT SCENE');
            else
                skipNextScene = 0;
                fprintf('PRE-TRIAL FIXATION: mean distance from center was %4f', fixShiftList(preFixCount,1));
            end
            % raw head - last 250 samples - this is only useful for the
            % contingent gaze code. pre-contingent gaze the fixations are a much better way 
            if length(pitch)>250
                meanHeadPitch = mean( pitch(end-250:end) ) - 90;
                meanHeadYaw = mean( yaw(end-250:end) ) - 180;
                meanHeadPitchList(preFixCount,1) = meanHeadPitch;
                meanHeadYawList(preFixCount,1) = meanHeadYaw;
                [headShiftList(preFixCount,1), foo] = distance(meanHeadYaw, meanHeadPitch, params.avgPreTrialX, params.avgPreTrialY); 
            else
                meanHeadPitch = nan;
                meanHeadYaw = nan;
                meanHeadPitchList(preFixCount,1) = meanHeadPitch;
                meanHeadYawList(preFixCount,1) = meanHeadYaw;
                headShiftList(preFixCount,1) = nan;
            end
            
            % raw gaze - last 250 samples
            if length(gaze_pitch_sphere)>250
                meanGazePitch = mean( gaze_pitch_sphere(end-250:end) ) - 90;
                meanGazeYaw = mean( gaze_yaw_sphere(end-250:end) ) - 180;
                meanGazePitchList(preFixCount,1) = meanGazePitch;
                meanGazeYawList(preFixCount,1) = meanGazeYaw;
                [gazeShiftList(preFixCount,1), foo] = distance(meanGazeYaw, meanGazePitch, params.avgPreTrialX, params.avgPreTrialY); 
            else
                meanGazePitch = nan;
                meanGazeYaw = nan;
                meanGazePitchList(preFixCount,1) = meanGazePitch;
                meanGazeYawList(preFixCount,1) = meanGazeYaw;
                gazeShiftList(preFixCount,1) = nan;
            end
            
            
            if params.driftCorrection == 1 % if doing drift correction, set the amount to shift gaze points
%                 preshiftX = meanXshift;
%                 preshiftY = meanYshift;

                preshiftX = fixShiftX;
                preshiftY = fixShiftY;
                
%                 preshiftX = viewportNorm2dva(fixShiftX);
%                 preshiftY = viewportNorm2dva(fixShiftY);
                
            end
        end
    end
    
    %set current scene image file, either .jpeg or .jpg
    if exist([paths.projectStimDir currentScene '.jpeg'])
        imFile = [paths.projectStimDir currentScene '.jpeg']; 
    elseif exist([paths.projectStimDir currentScene '.jpg'])
        imFile = [paths.projectStimDir currentScene '.jpg'];
    else
        imFile = '';
    end
    if exist(imFile) % if the image file exists plot things on it
        myIm = imread(imFile); % load image as matrix
        resizeIM = imresize(myIm,[params.imDimY,params.imDimX]);
        
        %Plot fixation points with spread
        if params.plotHeadFixationInfo == 1
            try
                [point_spread] = plotFix(resizeIM, mean_fix_yaw, mean_fix_pitch, rawSceneData.rawHeadYaw,rawSceneData.rawHeadPitch, begin_fix, length_fix-1, subjectName, currentSceneText, fileID, fix_duration, paths, params);
            catch 
            end
        end
        
        if params.plotCombinedGazeFlag == 1
            plotCombinedGaze(resizeIM, yaw, pitch, subjectName, currentSceneText, confidence, paths, gaze_yaw_sphere, gaze_pitch_sphere, mean_fix_yaw, mean_fix_pitch, begin_fix, length_fix, fileID, params)
        end
        if params.plotHeadRawFlag == 1 % Plot raw HMD Direction center
            plotHeadRaw(resizeIM, yaw, pitch, subjectName, currentSceneText, confidence, paths, params)
        end
        %Plot raw gaze points with spread
        if params.plotGazeRawFlag==1 %Plot Raw Fixation points
            plotGazeRaw(resizeIM, gaze_yaw_sphere, gaze_pitch_sphere, subjectName, currentSceneText, confidence, paths, params)
        end
        %Plot fixation points with spread
        if params.plotFixFlag == 1
            try
                [point_spread] = plotFix(resizeIM, mean_fix_yaw, mean_fix_pitch, gaze_yaw_sphere, gaze_pitch_sphere, begin_fix, length_fix, subjectName, currentSceneText, fileID, fix_duration, paths, params);
            catch 
            end
        end
    else % skip plotting if no image file
        fprintf('Skipping Scene: %s\n (no image file)', currentScene);%log current scene to console
    end
    metaData.(['scene' currentSceneText]).included = 1;
    %clear variables before next loop
    clearvars -except junkSceneList percentScannedList lowConfidenceCount sceneCount sceneProcessedList timeList fpsList uniqueID metaData avgFixDist preFixCount meanRawConfidenceList xShiftList yShiftList fixShiftList excludeScenes ...
        excludedScenesPar dataPath subjectName dataImport meanHeadPitch meanHeadPitchList meanHeadYawList meanHeadYaw headShiftList...
        parData sceneNames sceneChangeList fileID timeList paths params subjectName parIndex gazeShiftList meanGazeYawList meanGazePitchList skipNextScene...
        preshiftX preshiftY excludedPretrialCount minSampleCount eyeTrackerFailCount excludedScannedCount;
    
    close all;
end
meanTimeList = mean(timeList);
metaData.excludedPretrialCount = excludedPretrialCount;
metaData.minSampleCount = minSampleCount;
metaData.eyeTrackerFailCount = eyeTrackerFailCount;
metaData.percentScannedList = percentScannedList;
metaData.excludedScannedCount = excludedScannedCount;
metaData.lowConfidenceCount = lowConfidenceCount;
metaData.meanRawConfidenceList = meanRawConfidenceList;
metaData.meanTimeList = meanTimeList;
metaData.timeList = timeList;
metaData.fpsList = fpsList;
meanFPS = mean(fpsList);
metaData.meanFPS = meanFPS;

% if params.avgPreTrial == 1
%     metaData.xShiftList = xShiftList;
%     metaData.yShiftList = yShiftList;
%     metaData.fixShiftList = fixShiftList;
%     metaData.meanXShift = mean(xShiftList);
%     metaData.meanYShift = mean(yShiftList);
%     metaData.meanFixShift = mean(fixShiftList);
%     metaData.meanHeadPitch = mean(meanHeadPitch);
%     metaData.meanHeadYaw = mean(meanHeadYaw);
%     metaData.headShiftList = headShiftList;
%     metaData.meanHeadPitchList = meanHeadPitchList;
%     metaData.meanHeadYawList = meanHeadYawList;
%     %new sample based
%     metaData.meanGazePitchList = meanGazePitchList;
%     metaData.meanGazeYawList = meanGazeYawList;
%     metaData.gazeShiftList = gazeShiftList;
% end

metaData.sceneProcessedList = sceneProcessedList;
fprintf(fileID, '\n\nMean Time in Scene: %.2f seconds\n', meanTimeList);
fprintf(fileID, '\n\nMean FPS: %.2f frames/second\n', meanFPS);

if params.avgPreTrial == 1
    try
        fprintf(fileID, '\n\nAvg Pre-Trial Error: %3f degrees\n', mean(fixShiftList));
    catch
        
    end
end
fprintf(fileID, '\n\nEnd Time: %s\n', datetime);
save( [paths.projectMetaDataDir subjectName] , 'metaData' );
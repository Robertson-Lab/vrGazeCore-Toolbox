function [time,gazeX_deg,gazeY_deg,gaze_ecc,confidence,gaze_pitch_sphere,gaze_yaw_sphere,validIdxs,removedIdxs]=interpolSamples(gaze_yaw_sphere,gaze_pitch_sphere,time, gazeX_deg, gazeY_deg, gaze_ecc, confidence, parDataScene, rawSceneData, paths,params)
    % DEPRICATED function: not actively maintained or tested & not available in Python version
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

end
function [mean_fix_yaw,mean_fix_pitch,fix_duration,length_fix,start_time, end_time,begin_fix_idx,end_fix_idx,mean_fix_gaze_ecc, yaw_sphere, pitch_sphere] = calculateFixations(removedIdxs, gazeX_deg,gazeY_deg,yaw_sphere,pitch_sphere,time,gaze_ecc,fileID,paths,subjectName,currentSceneText,params,rawData)
%Calculates fixations based on mean absolute deviation
%This is the meat of the findFixations pipeline
% Inputs:
%   - removedIdxs = indexs removed due to a filter are marked and should be
%   considered when calculating fixations
%   - gazeX_deg = eye yaw in degrees (viewport space)
%   - gazeY_deg = eye pitch in degrees (viewport space)
%   - yaw_sphere = yaw of fixType in question in sphere space
%   - pitch_sphere =  pitch of fixType in question in sphere space
%   - time = experiment clock of current trial
%   - gaze_ecc = distance (degs) of eye from viewport center in a given scene
%   - fileID = which number scene are we on
%   - paths = all file paths
%   - subjectName = current participant ID
%   - currentSceneText = scene name + the order of presentation (e.g. the 3rd scene would be scenename003)
%   - params = all current params set by coreParams
%   - rawData = raw output from Unity
% 
% Outputs:
%   - mean_fix_yaw = matrix of mean yaw for each fixation point
%   - mean_fix_pitch = matrix of mean pitch for each fixation point
%   - fix_duration = time that a fixation lasted (secs)
%   - length_fix = length of a fixation in # of data points (row idxs)
%   - start_time = start time of a given fixation on the experiment clock
%   - end_time = end time of a given fixation on the experiment clock
%   - begin_fix_idx = start idx of a fixation
%   - end_fix_idx = end idx of a fixation
%   - mean_fix_gaze_ecc = matrix of mean gaze ecc for each fixation (how from viewport center was the fix point)

% tlb adding these - will probably want to clean this
%   - yaw_sphere = yaw with interpolated samples replaced by  gaze centroid
%   - pitch_sphere = pitch with interpolated samples replaced by gaze centroid


% Also exports fixations as .mat file for the current subject
% generalize this function so that it can be used for gaze (still), head,
% eye
%
% TLB 4-25-19 changing variable names, commenting, and QC-ing script
%  JSM note to add!!! - write out another column that is fixation Index

%% Distance Calculation (orthodromic)
%computes the lengths, arclen, of the great circle arcs connecting pairs of points on the surface of a sphere.

% move back to lat lon space for distance fn
if params.fixType == 1 % raw HMD yaw/pitch data is already in lat long space
    yaw_sphere = yaw_sphere-180; % gaze yaw in -180 to 180 space
    pitch_sphere = pitch_sphere-90; % gaze pitch in -90 to 90 space
end

for i = 1:length(yaw_sphere)-1
    [dis(i,1),az(i,1)] = distance(yaw_sphere(i), pitch_sphere(i), yaw_sphere(i+1), pitch_sphere(i+1)); %find orthodromic distance from one gaze point to the next
end

%% Velocity Calculation
for i = 1:length(time)-1
    deltaTime(i,1) = time(i+1) - time(i);
end
vel = dis./deltaTime; %find velocity of a given gaze point

%% Acceleration Calculation
% tlb - maybe useful for something?? idk
for i = 1:length(vel)-1
    deltavel(i,1) = vel(i+1) - vel(i);
end

accTime = deltaTime(1:end-1);
acc = deltavel./accTime;

%% MAD Calculation  
%compute Mean Absolute Deviation of distance 
%Peterson(2016): "Windows with a MAD less than 50�/s were classified as
%potential fixations, with consecutive qualifying windows concatenated into longer potential fixations"; 

swTime = time(2:end); %time adjusted for the sliding window --> removing the sample that gets removed in a distance calculation
swSize = 0.1; %size of the sliding window in secs

%determine start and end idxs based on the amount of time we can test
startSWIdx = find(abs(swTime - (swTime(1) + swSize/2)) == min(abs(swTime - (swTime(1) + swSize/2)))); 
endSWIdx = find(abs(swTime - (swTime(end) - swSize/2)) == min(abs(swTime - (swTime(end) - swSize/2))));
madVel = [];

for i = startSWIdx:endSWIdx

    %creates a variable width mad window based on time
    %find the point approx 0.1ms around current point --> 50ms ahead and 50ms behind
    beforeIdx = find(abs(swTime - (swTime(i) - swSize/2)) == min(abs(swTime - (swTime(i) - swSize/2))));
    afterIdx = find(abs(swTime - (swTime(i) + swSize/2)) == min(abs(swTime - (swTime(i) + swSize/2))));
    madVel(i,1) = mad(vel(beforeIdx:afterIdx));

end

%% Calculate the Fixations!
%%%%% Windows with a MAD less than x deg/s were classified as potential fixations,jsm - I got help from here: https://www.mathworks.com/matlabcentral/answers/224769-how-to-get-the-starting-and-ending-index-of-repeated-numbers-in-an-array
% Identify potential fixations and saccades

% ARM - figure out how to add fixSpatialDist and fixTempDist to coreParams. figure out how to reference those variables here
switch params.fixType
    case 1 % GAZE FIXATIONS - if we want to calculate gaze fixation use madVel
        fixCalcData = madVel(startSWIdx+1:end); %remove points not included in the mad sliding window
        fixSpatialDist = 2; % if fixations are within this # of deg, group together
        fixTempDist = 0.15; % if fixations are within fixTimeLimit of each other group togethewr
    case 2 % HEAD FIXATIONS - if we want to calculate head fixation, use velocity
        fixCalcData = vel;
        fixSpatialDist = 2;
        fixTempDist = 0.10;
end

potential_fixations = find(fixCalcData<params.minMad);%find a list of MAD idxs of potential fixations

%using the idxs of potential_fixations, find the spacing between idxs by taking the difference
potential_fix_diff=diff(potential_fixations);%equal to 1 when consecutive, >1 when not consecutive (saccade possibly)

% find all non-consecutive potential fixes, will be potential_fix_diff will be > 1 when changing b/w two discrete gaze points
potential_fix_start= [1 ; (potential_fix_diff > 1)]; % first point is always a fix start, add 1 to beginning

for i = 1:length(potential_fix_start)-1 % tlb building a custom difference function to find the end fixation point
    % if the next point is not neighboring, the point will have a value one
    if (potential_fix_start(i+1) > 0 || potential_fix_start(i) > 0)
         %a fix start followed by a neighboring pt will have a value of -1
         %a fix start followed by another fix start will have a value of 0
         %a fix start PRECEDED by a neighboring pt will have a value of 1
        potential_fix_end(i,1) = (potential_fix_start(i+1) - potential_fix_start(i) > -1);
    else
        potential_fix_end(i,1) = 0;
    end
end

%potential_fix_end = diff([1 ; potential_fix_diff]) > 0; % first point can never truly be the end of a fixation
begin_fix_idx= potential_fixations(find(potential_fix_start.'));%find beginning of each group, 'saccades', add 1 at beginning b/c first point is always beginning a fixation
end_fix_idx= potential_fixations(find([potential_fix_end.',true]));%find end of each group
length_fix=1+end_fix_idx-begin_fix_idx;%find length of each group

if params.fixType == 1

    time = time(startSWIdx+1:endSWIdx); % 1 point removed from beginning during distance calc 3 points removed from beginning and end of mad calculation

else % only 1 point is removed from original data when calculating vel
    begin_fix_idx = begin_fix_idx + 1; 
    end_fix_idx = end_fix_idx + 1;
end


%Calculate Centroid (center of mass) of points ON A SPHERE - LAT LON - from each fixation

removedIdxs = removedIdxs(2:end); % adjust the matrix containing which samples are valid to match the distance calculation 

for i=1:length(begin_fix_idx)
    
    %tlb adding this because i don't think we truly want to consider interpolated samples for calculating fixations

    currentPitchList = pitch_sphere(begin_fix_idx(i):end_fix_idx(i)); % yaw points from start of fix event to end of fix event
    currentYawList = yaw_sphere(begin_fix_idx(i):end_fix_idx(i)); % pitch points from start of fix event to end of fix event

    if params.unityProjectVersion == 1
        try 
            validIdxs = find(removedIdxs(begin_fix_idx(i):end_fix_idx(i))==0); % find idxs of current valid samples
            invalidIdxs = find(removedIdxs(begin_fix_idx(i):end_fix_idx(i))==1);

            [mean_fix_pitch(i),mean_fix_yaw(i)] = meanm(currentPitchList(validIdxs),currentYawList(validIdxs)); % mean pitch & yaw of the fix event

            if ~isempty(invalidIdxs) % if there are idxs that were from interpolation
                invalidIdxs = invalidIdxs + begin_fix_idx(i) - 1; % shift the idxs into alignment with the gaze coords
            end

            pitch_sphere(invalidIdxs) = mean_fix_pitch(i); % then assign the interpolated samples the mean gaze centroid
            yaw_sphere(invalidIdxs) = mean_fix_yaw(i);
        catch
            [mean_fix_pitch(i),mean_fix_yaw(i)] = meanm(currentPitchList,currentYawList); % mean pitch & yaw of the fix event
        end
    else
        
        [mean_fix_pitch(i),mean_fix_yaw(i)] = meanm(currentPitchList,currentYawList); % mean pitch & yaw of the fix event

    end

    start_time(i) = time(begin_fix_idx(i)); % time at which a fix event began
    end_time(i) = time(end_fix_idx(i)); % time of the fixation ended
    fix_duration(i) = end_time(i) - start_time(i);% amount of time between start of fixation and start of next non-fix event
end

potential_fix_count = length(begin_fix_idx);

%% Concatenate Fixations
%If fixations are very close temporally and spatially, add them together.
% Peterson (2016): Potential fixations separated by fewer than nine consecutive 
% invalid samples (150 ms) were concatenated if they were displaced by less than 1 deg
% with invalid samples assigned the mean position of the preceding potential fixation

%tlb 8-29-19 adding in method to account for invalid samples --> assigning
%mean position of preceding sample

concatenate = [];
concount = 0;
removeCount = 1;
for i=1:length(begin_fix_idx)-1
    
%     % find amount of time from invalid samples between the end of current fixation and
%     % beginning of next fixation
%     
%     %check to see if curr removed idx is between two potential fixation
%     %idxs and the two values are separated by less than 150ms, potential for concatenation
%     curr_removal_idxs = find(removedIdxs(:,1) >= end_fix_idx(i) & removedIdxs(:,1) < begin_fix_idx(i+1)); 
%     %removedIdxs(removeCount,1) >= end_fix_idx(i) && removedIdxs(removeCount,1) < begin_fix_idx(i+1) &&
    
    [dist_degrees,d_az] = distance(mean_fix_yaw(i),mean_fix_pitch(i), mean_fix_yaw(i+1), mean_fix_pitch(i+1)); %find distance between fixation i and i+1 (arclen degrees)
    if start_time(i+1)-end_time(i) < fixTempDist && dist_degrees < fixSpatialDist %% should base this on our accuracy?
        concatenate(i) = 1;
    else
        concatenate(i) = 0;
    end
end

i = 1;
while i<=length(concatenate) %consider the last point!
    if concatenate(i) == 1
        [mean_fix_pitch(i),mean_fix_yaw(i)] = meanm([mean_fix_pitch(i); mean_fix_pitch(i+1)],[mean_fix_yaw(i); mean_fix_yaw(i+1)]); %find the mean of the two pitch/yaw values
        mean_fix_yaw(i+1) = []; % remove next yaw value b/c has been averaged in
        mean_fix_pitch(i+1) = []; % remove next pitch value b/c has been averaged in
        fix_duration(i) = fix_duration(i) + fix_duration(i+1);%add durations, delete next one
        fix_duration(i+1) = []; %add durations
        length_fix(i+1) = []; % remove next duration b/c has been added
        start_time(i+1) = []; % remove next length fixation
        end_time(i) = []; % remove current end time
        begin_fix_idx(i+1) = []; % remove next fix idx b/c has been grouped
        end_fix_idx(i) = []; % remove current end idx
        concatenate(i)=[]; %remove current entry from concatenate
        concount = concount + 1;
    else
        i = i + 1; %otherwise look for next 1
    end

end
% log concatenation results.
if params.fixType == 1 % only want to do this for eye fix? tlb solve later
    fprintf(fileID,'GAZE CONCATENATION: concatenated %d out of %d potential fixations\n', concount, potential_fix_count);
    fprintf('GAZE CONCATENATION: concatenated %d out of %d potential fixations\n', concount, potential_fix_count);
end

%% Filter out first n and last n fixations
excludeAmount = params.excludeFirstNFix + 1; % need to add 1 to get proper idx

if length(begin_fix_idx) > params.excludeFirstNFix && params.fixType == 1 % make sure # of fixations is > # of wanted exclusions
    if params.excludeFirstNFix ~=0 && ~contains(currentSceneText, '_sanityTarget') % don't exclude fixations of the sanity scene --> tlb check this
        mean_fix_yaw = mean_fix_yaw(excludeAmount : length(mean_fix_yaw)-params.excludeLastNFix);
        mean_fix_pitch = mean_fix_pitch(excludeAmount : length(mean_fix_pitch)-params.excludeLastNFix);
        fix_duration = fix_duration(excludeAmount : length(fix_duration)-params.excludeLastNFix);
        length_fix = length_fix(excludeAmount : length(length_fix)-params.excludeLastNFix);
        start_time = start_time(excludeAmount : length(start_time)-params.excludeLastNFix);
        end_time = end_time((excludeAmount : length(end_time)-params.excludeLastNFix));
        begin_fix_idx = begin_fix_idx(excludeAmount : length(begin_fix_idx)-params.excludeLastNFix);
        end_fix_idx = end_fix_idx((excludeAmount : length(end_fix_idx)-params.excludeLastNFix));
    end
end

%% Validate fixations

% validity notes
%%% Potential fixations were labeled as valid fixations IF
%%% they were immediately preceded and followed by a likely saccade event
%%% (MAD of three preceding/succeeding samples greater than 100 degrees/s) AND
%%% displaced from the mean of the preceding and succeeding potential fixations by at least 1 degree.

% % %tlb 8-29-19 - then exclude same fixation from verification

% tlb 8-29-19 shift back into vel space by adding in removed pts from mad
% window
% 
% shift vel into mad space?

if params.fixValidation == 1
    prefix_pt_validation = begin_fix_idx; %- (s_window/2) - 1;
    postfix_pt_validation = end_fix_idx;%- (s_window/2) - 1;
    params.valWindow = 2; %number of points before a fixation to check

    % time(prefix_pt_validation(3)-10:prefix_pt_validation(3))
    % dis(prefix_pt_validation(3)-10:prefix_pt_validation(3))

    % 
    for i = 2:length(begin_fix_idx)
        try
            % going to assume this means using the mad numbers already given
            prefix_pts(i,:) =  fixCalcData(prefix_pt_validation(i)-val_window:prefix_pt_validation(i)-1);
            postfix_pts(i,:) = fixCalcData(postfix_pt_validation(i)+1:postfix_pt_validation(i)+val_window);

    %         prefix_pts(i) = mad(vel(prefix_pt_validation(i)-val_window:prefix_pt_validation(i)));
    %         postfix_pts(i) = mad(vel(postfix_pt_validation(i):postfix_pt_validation(i)+val_window));    
    % 
    %         
    %         postfix_pts(i) = mad(vel(postfix_pt_validation(4):postfix_pt_validation(4)+val_window));    
    %         prefix_pts(i) = mad(vel(prefix_pt_validation(4)-val_window:prefix_pt_validation(4)));
        catch

        end
    %     prefix_pts = mad_distance(prefix_pt_validation(i)-val_window:prefix_pt_validation(i));
    %     postfix_pts = mad_distance(postfix_pt_validation(i):postfix_pt_validation(i)+val_window);
    end


    try
        prefix_pts = prefix_pts(excludeAmount:end,:);
        postfix_pts = postfix_pts(excludeAmount:end,:);

        verifyLim = params.maxMad;


        validFix = intersect(find(all(prefix_pts(:,:)> verifyLim, 2)),find(all(postfix_pts(:,:) > verifyLim, 2)));

        %if you want the first fixation included, add it on
        validFix = [1 ; validFix];

        %then remove fixations that wouldn't qualify under MAD parameter
        mean_fix_yaw = mean_fix_yaw(validFix);
        mean_fix_pitch = mean_fix_pitch(validFix);
        fix_duration = fix_duration(validFix);
        length_fix = length_fix(validFix);
        start_time = start_time(validFix);
        end_time = end_time(validFix);
        begin_fix_idx = begin_fix_idx(validFix);
        end_fix_idx = end_fix_idx(validFix);
    catch

    end
end

%% Then filter out any fixations less than 100ms of time
if params.excludeFixDursLessThan ~= 0
    shortFixIdx = find(fix_duration < params.excludeFixDursLessThan);
    mean_fix_yaw(shortFixIdx) = [];
    mean_fix_pitch(shortFixIdx) = [];
    fix_duration(shortFixIdx) = [];
    length_fix(shortFixIdx) = [];
    start_time(shortFixIdx) = [];
    begin_fix_idx(shortFixIdx) = [];
    end_fix_idx(shortFixIdx) = [];
    end_time(shortFixIdx) = [];
end

%% Calculate Mean Fixation Eccentricity %JSM CHECK THIS

if params.fixType == 1 && length(begin_fix_idx > 1)
%     for i=1:length(begin_fix_idx)
%         mean_fix_gaze_ecc(i) = mean(gaze_ecc(begin_fix_idx(i):end_fix_idx(i)));
%     end
    mean_fix_gaze_ecc = []; % initialize matrix to prevent output error 
else
    mean_fix_gaze_ecc = []; % initialize matrix to prevent output error 
end


if isnan(fix_duration)
    h
end

% TLB - checking common fixation stats
fprintf('\nNum Fixations: %d\n', length(fix_duration));
fprintf('Avg Fixation Duration: %0.3f\n', mean(fix_duration));
fprintf('STD Fixation Duration: %0.3f\n', std(fix_duration));
fprintf('Min Fixation Duration: %0.3f\n', min(fix_duration));
fprintf('Max Fixation Duration: %0.3f\n\n', max(fix_duration));

%% Export Fixations
%Gather spherical median fixations, bring back to positive domain
if (params.fixType == 1 || params.fixType == 2) %only export this data if we are doing gaze
    %TLB - maybe change later if we want to integrate headFixations into findFixations
    
    mean_fix_yaw = mean_fix_yaw+180;
    mean_fix_pitch = mean_fix_pitch+90;
    xSphereMedian = mean_fix_yaw;
    ySphereMedian = mean_fix_pitch;
    startMedian = start_time;
    durationMedian = fix_duration;

    %Gather spherical raw gaze that corresponds to fixation times (spherical)
    yaw_sphere = yaw_sphere+180;
    pitch_sphere = pitch_sphere+90;
    xSphereRaw = [];
    ySphereRaw = [];
    startRaw = [];
    durationRaw = [];

    %Gather HMD raw gaze (degrees ~FoV) that corresponds to fixation times
    xHMDRaw = [];
    yHMDRaw = [];
%     for ifix = 1:length(begin_fix_idx)
%         for wfix = 1:length_fix(ifix)
%             xSphereRaw = [xSphereRaw yaw_sphere(begin_fix_idx(ifix)+wfix)];
%             ySphereRaw = [ySphereRaw pitch_sphere(begin_fix_idx(ifix)+wfix)];
%             startRaw = [startRaw time(begin_fix_idx(ifix)+wfix)];
%             durationRaw = [durationRaw (time(begin_fix_idx(ifix)+wfix + 1)) - time(begin_fix_idx(ifix)+wfix) ];
%             xHMDRaw = [xHMDRaw gazeX_deg(begin_fix_idx(ifix)+wfix)];
%             yHMDRaw = [yHMDRaw gazeY_deg(begin_fix_idx(ifix)+wfix)];
%         end
%     end

    %JSM CHECK THIS - put hardcoded vals in parms file (including time window above)
    %Gather Equirectangular median and raw fixations  %%%%%%%%%%%%%%%%%%%%
    [xEquirectMedian,yEquirectMedian] = wrapPointsEquirect( mean_fix_yaw, mean_fix_pitch, 360, 180, 0); 
    [xEquirectMedian, yEquirectMedian] = degreesToPixels(xEquirectMedian,yEquirectMedian,params.imDimX,params.imDimY);
    [xSphereRaw,ySphereRaw] = wrapPointsEquirect( xSphereRaw, ySphereRaw, 360, 180, 0); 
    [xEquirectRaw, yEquirectRaw] = degreesToPixels(xSphereRaw,ySphereRaw,params.imDimX,params.imDimY);
    
    % if ~exist( [paths.projectFixDataDir  subjectName '/'] )
    %     mkdir( [paths.projectFixDataDir  subjectName '/'] );
    % end

    %%Create structures  to save in .mat file

    %Fields for fixation data (in spherical, viewport, or equirectangular space)
    fixData = struct;
    fixData.gazeFixXSphere = xSphereMedian;
    fixData.gazeFixYSphere = ySphereMedian;

    fixData.gazeFixXEqui = xEquirectMedian;
    fixData.gazeFixYEqui = yEquirectMedian;

    fixData.gazeFixStart = start_time;
    fixData.gazeFixEnd = end_time;
    fixData.gazeFixDur = fix_duration;

    %Fields for raw gazepoints that underlie each fixation
    rawFixData = struct;
    
    rawFixData.filteredVPSGazeX = gazeX_deg;
    rawFixData.filteredVPSGazeY = gazeY_deg;
    rawFixData.filteredSSgazeX = yaw_sphere;
    rawFixData.filteredSSgazeY = pitch_sphere;
    
    rawFixData.gazeUnderlyingFixationX = xSphereRaw;
    rawFixData.gazeUnderlyingFixationY = ySphereRaw;

    rawFixData.viewportUnderlyingFixationX = xHMDRaw;
    rawFixData.viewportUnderlyingFixationY = yHMDRaw;

    rawFixData.gazeUnderlyingFixationStart = startRaw;
    rawFixData.gazeUnderlyingFixationDuration = durationRaw;

    save([paths.subjectFixMatDir currentSceneText ],'rawData', 'fixData','rawFixData');

end
end
function [sacc_durs,length_sacc,sacc_start_times, sacc_end_times,begin_sacc_idx,end_sacc_idx,mean_sacc_gaze_ecc] = calculateSaccades(xGazeVPS,yGazeVPS,time,gaze_ecc,fileID,paths,subjectName,currentSceneText,params,rawData)
%Calculates fixations based on mean absolute deviation
%This is the meat of the findFixations pipeline
% Inputs:
%   - xGazeVPS = eye yaw in degrees (viewport space)
%   - yGazeVPS = eye pitch in degrees (viewport space)
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
% TLB need to figure out mean_fix_yaw and stuff 
%   - mean_fix_yaw = matrix of mean yaw for each fixation point
%   - mean_fix_pitch = matrix of mean pitch for each fixation point
%   - sacc_durs = time that a fixation lasted (secs)
%   - length_sacc = length of a saccade in # of data points (row idxs)
%   - start_time = start time of a given fixation on the experiment clock
%   - end_time = end time of a given fixation on the experiment clock
%   - begin_fix_idx = start idx of a fixation
%   - end_fix_idx = end idx of a fixation
%   - mean_fix_gaze_ecc = matrix of mean gaze ecc for each fixation (how from viewport center was the fix point)
%
% Also exports fixations as .mat file for the current subject
% generalize this function so that it can be used for gaze (still), head,
% eye
%
% TLB 4-25-19 changing variable names, commenting, and QC-ing script
%  JSM note to add!!! - write out another column that is fixation Index

%% Distance Calculation (orthodromic)
%computes the lengths, arclen, of the great circle arcs connecting pairs of points on the surface of a sphere.

for i = 1:length(xGazeVPS)-1
    [sacc_dis(i)] = pdist2([xGazeVPS(i) xGazeVPS(i+1)],[yGazeVPS(i) yGazeVPS(i+1)]);
end

%% Velocity Calculation
for i = 1:length(time)-1
    deltaTime(i,1) = time(i+1) - time(i);
end
sacc_vel = sacc_dis'./deltaTime;

%% MAD Calculation
s_window = 6; %set sliding window size in samples (even number). Value 6 will give a 7 sample window, 3 before and 3 after each point. 
for i = 1+(s_window/2):length(sacc_vel)-s_window
    sacc_mad_distance(i,1) = mad(sacc_vel((i-(s_window/2)):(i+(s_window/2))));
end

%% Acceleration Calculation
% tlb - maybe useful for something?? idk
for i = 1:length(sacc_vel)-1
    delta_sacc_vel(i,1) = sacc_vel(i+1) - sacc_vel(i);
end
sacc_acc = delta_sacc_vel./deltaTime(1:end-1);

%% Find saccades
%%% Potential fixations were labeled as valid fixations IF
%%% they were immediately preceded and followed by a likely saccade event
%%% (MAD of three preceding/succeeding samples greater than 100 degrees/s) AND
%%% displaced from the mean of the preceding and succeeding potential fixations by at least 1 degree. 

% Peterson (2016): Saccades were classified as events where eye velocity was 
% greater than 22deg/s and eye acceleration exceeded 4000deg/s^2 

% TLB going to make both methods and see what works

% mad > 100m/s acceleration has to be > 4000m/s^2 and velocity > 22m/s at the same
% point in time to be considered a saccade
%  potential_saccades = intersect(find(sacc_acc > 4000),find(sacc_vel(1:end-1)>22));

potential_saccades = find(sacc_mad_distance > 100);
potential_sacc_diff= diff(potential_saccades);%equal to 1 when consecutive, >1

% find all non-consecutive potential fixes, will be potential_fix_diff will be > 1 when changing b/w two discrete gaze points
potential_sacc_start= [1 ; (potential_sacc_diff > 1)]; % first point is always a fix start, add 1 to beginning
%any 0 will be a change


%begin_sacc_idx = begin_sacc_idx(acc(begin_sacc_idx) >= 4000); %filter based on acceleration

for i = 1:length(potential_sacc_start)-1 %building a custom difference function to find the end fixation point
    % if the next point is not neighboring, the point will have a value one
    if (potential_sacc_start(i+1) > 0 || potential_sacc_start(i) > 0)
         %a fix start followed by a neighboring pt will have a value of -1
         %a fix start followed by another fix start will have a value of 0
         %a fix start PRECEDED by a neighboring pt will have a value of 1
        potential_sacc_end(i,1) = (potential_sacc_start(i+1) - potential_sacc_start(i) > -1);
    else
        potential_sacc_end(i,1) = 0;
    end
end

begin_sacc_idx= potential_saccades(find(potential_sacc_start.'));%find beginning of each group, 'saccades'
end_sacc_idx= potential_saccades(find([potential_sacc_end.',true]));%find end of each group
length_sacc = 1 + end_sacc_idx - begin_sacc_idx;

begin_sacc_idx = begin_sacc_idx(find(length_sacc ~= 0));
end_sacc_idx = end_sacc_idx(find(length_sacc ~= 0));

sacc_gaze_start_idx = begin_sacc_idx+ 1; 
sacc_gaze_end_idx = end_sacc_idx + 1;

for i = 1:length(begin_sacc_idx)
    %sacc_mean_mad(i) = mean(mad_distance(begin_sacc_idx(i):end_sacc_idx(i))); % find the mad of each point in the saccade
    sacc_start_times(i) = time(sacc_gaze_start_idx(i));
    sacc_end_times(i) = time(sacc_gaze_end_idx(i));
    sacc_durs(i) = time(sacc_gaze_end_idx(i)) - time(sacc_gaze_start_idx(i)); % find the mean time of the saccade
end


%% Do we want to exclude saccades??
% excludeAmount = params.excludeFirstNFix + 1; % need to add 1 to get proper idx
% if length(begin_fix_idx) > params.excludeFirstNFix && params.fixType == 1 % make sure # of fixations is > # of wanted exclusions
%     if params.excludeFirstNFix ~=0
%         mean_fix_yaw = mean_fix_yaw(excludeAmount : length(mean_fix_yaw)-params.excludeLastNFix);
%         mean_fix_pitch = mean_fix_pitch(excludeAmount : length(mean_fix_pitch)-params.excludeLastNFix);
%         fix_duration = fix_duration(excludeAmount : length(fix_duration)-params.excludeLastNFix);
%         length_fix = length_fix(excludeAmount : length(length_fix)-params.excludeLastNFix);
%         start_time = start_time(excludeAmount : length(start_time)-params.excludeLastNFix);
%         begin_fix_idx = begin_fix_idx(excludeAmount : length(begin_fix_idx)-params.excludeLastNFix);
%     end
% end

%% (filter out short saccades??)
% if params.excludeFixDursLessThan ~= 0
%     shortFixIdx = find(fix_duration < params.excludeFixDursLessThan);
%     mean_fix_yaw(shortFixIdx) = [];
%     mean_fix_pitch(shortFixIdx) = [];
%     fix_duration(shortFixIdx) = [];
%     length_fix(shortFixIdx) = [];
%     start_time(shortFixIdx) = [];
%     begin_fix_idx(shortFixIdx) = [];
%     end_fix_idx(shortFixIdx) = [];
%     end_time(shortFixIdx) = [];
% end

%% Calculate Mean Saccade Eccentricity %do we want this
for i=1:length(begin_sacc_idx)
    mean_sacc_gaze_ecc(i) = mean(gaze_ecc(begin_sacc_idx(i):end_sacc_idx(i)));
end

%% Export Fixations
%Gather spherical median fixations, bring back to positive domain

% if ~exist( [paths.projectSaccadeDataDir  subjectName '/'] )
%     mkdir( [paths.projectSaccadeDataDir  subjectName '/'] );
% end

%%Create structures  to save in .mat file
    
%Fields for saccade data (will make better moving along - TLB) 
saccData = struct;
% find the idxs of saccades in the raw data for aligning w/ head data
[~, saccData.rawSaccStartIdxs] = ismember(xGazeVPS(sacc_gaze_start_idx),rawData.rawViewportX_deg); 
[~, saccData.rawSaccEndIdxs] = ismember(xGazeVPS(sacc_gaze_end_idx), rawData.rawViewportX_deg); 
saccData.saccStartIdxs =  sacc_gaze_start_idx;
saccData.saccEndIdxs = sacc_gaze_end_idx;
saccData.saccStartTimes = sacc_start_times;
saccData.saccEndTimes = sacc_end_times;
saccData.saccDurations = sacc_durs;
saccData.meanSaccGazeEcc = mean_sacc_gaze_ecc;

save([paths.subjectSaccadeDir erase(currentSceneText,'.mat') ],'saccData');

end
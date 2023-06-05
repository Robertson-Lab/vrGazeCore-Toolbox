function [parData,sceneNames, sceneChangeList] = processRawData(subjectName, params, paths)

% processRawData is a function that loads a data file of a participant and
% organizes it into a matrix to be processed by findFixations
%
% Inputs:
%   - subjectName: the name of the current participant to process (string)
%   - params: parameters set by coreParams script
%   - paths: paths to the vrGazeCore directories
%
% Outputs:
%   - parData: a matrix containing the current participant's data to be
%   processed by find fixations
%   - sceneNames: scene name associated with each collected sample. used
%   to find the start of each scene (stored in sceneChangeList)
%   - sceneChangeList: a vector of indices indicating the sample at the start of
%   a given trial scene. calculated by the difference in sceneNames
%
% written by TLB 10/23/19

%% Load Scene Gaze Data
try
    dataImport = importdata([char(paths.projectRawDataDir) char(subjectName) '.txt'], ','); %import raw data (8-10 columns)
catch
    fprintf(['\nThe data file of ' char(subjectName) ' does not exist! Skipping to the next file.']);
    return;
end_try_catch


%% Process Gaze data

% make ds of nans with 7 columns. columns represent the following:
%   - col1: experiment clock time
%   - col2: HMD yaw
%   - col3: HMD pitch
%   - col4: HMD roll
%   - col5: eye-tracking X coordinates
%   - col6: eye-tracking Y coordinates
%   - col7: eye-tracking confidence

parData = nan(size(dataImport.data,1), 7); % start an empty mat of nans
parData(:,1:4) = dataImport.data(:,1:4); % load in exp time and HMD data

if params.headsetType == 3 % headtracking only (oculus go)
    parData(:,5:end) = [];
    
else %otherwise eyetracking
    
    if params.headsetType == 2 &&  params.gazeType == 1 %  Vive Eye 3D tracking 
        switch params.useEye
            case 0 %use right eye only, col 9 = right eye confidence
                dataCols = [12, 13, 9]; 
            case 1 % use left eye only, col 10 = left eye confidence
                dataCols = [14, 15, 10];
            case 2 % best eye
                dataCols = [12, 13, 14, 15, 9, 10]; % right eye X Y, left eye X Y, right eye conf, left eye conf 
            case 3 % avg two eyes
                dataCols = [12, 13, 14, 15, 9, 10]; % right eye X Y, left eye X Y, right eye conf, left eye conf 
        end  
    else % DK2 or Vive Pro (w/ pupil trackers)

        % set which columns to grab data from based on useEye
        switch params.useEye
            case 0 %use right eye only, col 9 = right eye confidence
                dataCols = [5, 6, 9]; 
            case 1 % use left eye only, col 10 = left eye confidence
                dataCols = [7, 8, 10];
            case 2
                dataCols = [5, 6, 7, 8, 9, 10]; % right eye X Y, left eye X Y, right eye conf, left eye conf 
            case 3
                dataCols = [5, 6, 7, 8, 9, 10];
        end

    end

    if size(dataImport.data,2) >= 10  %if the data has more than 10 rows it is binocular, else if 7 rows monocular

        if params.useEye == 2 % "Best Eye"

            parData(:,5) = dataImport.data(:,dataCols(1)); %set normX from eye 0
            parData(:,6) = dataImport.data(:,dataCols(2)); %set normY from eye 0  
            parData(:,7) = dataImport.data(:,dataCols(5)); % set eye0 confidence to 7th column as in monocular data format

            for i = 1:size(dataImport.data) %loop through every timepoint
                if dataImport.data(i,dataCols(6)) > dataImport.data(i,dataCols(5)) % if eye 1 conf > eye 0 then replace eye 0 with 1 norm pos and conf
                    parData(i,5) = dataImport.data(i,dataCols(3)); %set normX from eye 1
                    parData(i,6) = dataImport.data(i,dataCols(4)); %set normY from eye 1
                    parData(i,7) = dataImport.data(i,dataCols(6)); %set confidence from eye 1
                end
            end

        elseif params.useEye == 3 % Center Point 

            parData(:,7) = min( cat(2, dataImport.data(:,dataCols(5)), dataImport.data(:,dataCols(6))) ,[],2 ); % set confidence to the minimum of the 2 eyes
            parData(:,5) = ( dataImport.data(:,dataCols(1))+dataImport.data(:,dataCols(3)) ) / 2; %set normX from average of eye 0 and 1
            parData(:,6) = ( dataImport.data(:,dataCols(2))+dataImport.data(:,dataCols(4)) ) / 2; %set normY from average of eye 0 and 1

        else % otherwise use a single eye

            parData(:,5) = dataImport.data(:,dataCols(1)); %set normX from specified eye
            parData(:,6) = dataImport.data(:,dataCols(2)); %set normY from specified eye
            parData(:,7) = dataImport.data(:,dataCols(3)); % set eye confidence to 7th column as in monocular data format

        end  
    else

        parData = dataImport.data(:,1:7); % monocular pupil data

    end

    % Vive Eye 2D tracking needs to be normalized to same field as pupil
    if params.headsetType == 2 && params.gazeType == 0 
        parData(:,5:6) = (parData(:,5:6) + 1)/2; % bring positive, then half
    end
end

%% SPLIT DATA BY SCENE
sceneNames = dataImport.textdata(:,1);
if params.concatSanity == 1 % change all pre-trial fixation scenes to have the correct same name for plotting, yes we do want to do this actually don't remove jsm 4/22/19
    sceneNames(find(contains(sceneNames, '_sanityTarget360'))) = {'_sanityTarget360'};
end

for l=1:length(sceneNames)-1
    scenechange(l)=strcmp(sceneNames(l),sceneNames(l+1));
end

%%  Set how the scene change points are calculated based on the unity project version

   
    if find(scenechange==0) > 1
        sceneChangeList = [find(scenechange==0) length(sceneNames)];%find indices where a scene changes, add 1 to beginning since that is the first scene
        
    else
        sceneChangeList = [1 find(scenechange==0) length(sceneNames)];%find indices where a scene changes, add 1 to beginning since that is the first scene
    end

end


function [fixationStartTimes, fixationEquiX, fixationEquiY, fixationDurations] = getParticipantFixations(subjectName, currentScene, fileID, paths, params)

%         subjectName = char (inputSubjects(parIndex,:) ); % get subject name
    dirCurrentScene = dir([paths.projectFixMatDir '/' subjectName '/' currentScene '*.mat']);

    % for the current subject, return x,y coordinates on the equirectangular 
    % image along with durations of each fixation
    fixationEquiX = [];
    fixationEquiY = [];
    fixationDurations = [];

    if isempty(dirCurrentScene) % skip if this subject didn't see this scene(if mat file DNE)
        %skip
        fprintf(fileID,'\n- - - - - - - - -%s did not see scene, skipping', subjectName);%log current subject to console
        fprintf('\n- - - - - - - - -%s did not see scene, skipping', subjectName);%log current subject to console
        
        fixationStartTimes = [];
        fixationEquiX = [];
        fixationEquiY = [];
        fixationDurations = [];
        return
    end
    
    % otherwise the file exists and we can load
    fprintf(fileID,'\n-+-+-+-+- Adding %s fixations to list', subjectName);%log current scene# to console
    fprintf('\n-+-+-+-+- Adding %s fixations to list', subjectName);%log current scene# to console

    loadFile = [paths.projectFixMatDir subjectName '/' dirCurrentScene(end).name] ;
    loadedfile = load(sprintf(loadFile));

    sceneStartTime = loadedfile.rawData.rawTime(1);
    fixationStartTimes = loadedfile.fixData.gazeFixStart;

    % set variables for output --> this could be cleaner as we're wasting
    % memory on holding two variables

    fixationStartTimes = fixationStartTimes - sceneStartTime; % this moves the fixation start times to be times relative to 0 (the scene start time)
    fixationEquiX = loadedfile.fixData.gazeFixXEqui;
    fixationEquiY = loadedfile.fixData.gazeFixYEqui;
    fixationDurations = loadedfile.fixData.gazeFixDur; %create variable to hold current subject's fixation durations

end
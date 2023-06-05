%% Script to take the subject list and plot all of the meta data from those subejcts


loadCoreParams;

dataStruct = struct;
%% load everyone into a struct
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    subjectName = char (params.subjectNames{parIndex,:} ); % get subject name
    cohortName = subjectName(1:end-3);
    dataStruct.sub{parIndex,1} = load([paths.projectDir 'data/meta/' subjectName '.mat']);
    dataStruct.subjectName{parIndex,1} = subjectName;
    dataStruct.cohortName{parIndex,1} = cohortName; 
end
    
%% Loop through the struc tand plot it all 


%%plot Confidence Scene Viewing
figure('DefaultAxesFontSize',13)
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = dataStruct.sub{parIndex,1}.metaData.meanRawConfidenceList;
    plot(  plotThis(plotThis~=0)    ,'LineWidth',2);
    clear plotThis
end

ylim([0 1])
ylabel('Pupil Confidence')
xlabel('Looking Around Sphere Scene #')
legend(char(dataStruct.subjectName'))
tidyFigs
title('Mean Raw Pupil Confidence');

%%plot x drift of eye
figure('DefaultAxesFontSize',13)
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = dataStruct.sub{parIndex,1}.metaData.xShiftList;
    plot(  plotThis(plotThis~=0)    ,'LineWidth',2);
    clear plotThis
end
ylim([-15 15])
ylabel('Drift - PreTrial - Degrees on Sphere')
xlabel('pretrial gaze scene #')
legend(char(dataStruct.subjectName'));
tidyFigs
title('Mean gaze drift (x dimension)');

%%plot y drift of eye
figure('DefaultAxesFontSize',13)
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = dataStruct.sub{parIndex,1}.metaData.yShiftList;
    plot(  plotThis(plotThis~=0)    ,'LineWidth',2)
    clear plotThis
end
ylim([-5 5])
ylabel('Drift - PreTrial - Degrees on Sphere')
xlabel('pretrial gaze scene #')
legend(char(dataStruct.subjectName'));
tidyFigs
title('Mean gaze drift (y dimension)');

%%plot total drift of Both Eyes
figure('DefaultAxesFontSize',13)
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = dataStruct.sub{parIndex,1}.metaData.fixShiftList;
    plot(   plotThis(plotThis~=0)   ,'LineWidth',2)
    clear plotThis
end
ylim([0 30])
ylabel('Drift - PreTrial - Degrees on Sphere')
xlabel('pretrial gaze scene #')
legend(char(dataStruct.subjectName'))
tidyFigs
title('Mean AVG Eye drift');





%%plot total head drift
figure('DefaultAxesFontSize',13)
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = dataStruct.sub{parIndex,1}.metaData.headShiftList;
    plot(   plotThis(plotThis~=0)   ,'LineWidth',2)
    clear plotThis
end
ylim([-5 5])
ylabel('Head Drift - PreTrial - Degrees on Sphere')
xlabel('pretrial gaze scene #')
legend(char(dataStruct.subjectName'))
tidyFigs
title('Mean AVG Head drift, last 250 samples of pre-trial scene');

%%plot total mean head yaw
figure('DefaultAxesFontSize',13)
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = dataStruct.sub{parIndex,1}.metaData.meanHeadYawList;
    plot(   plotThis(plotThis~=0)   ,'LineWidth',2)
    clear plotThis
end
ylim([-5 5])
ylabel('Head Drift - PreTrial - Degrees on Sphere')
xlabel('pretrial gaze scene #')
legend(char(dataStruct.subjectName'))
tidyFigs
title('Mean AVG Head Yaw, last 250 samples of pre-trial scene');




%%plot total mean head pitch
figure('DefaultAxesFontSize',13)
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = dataStruct.sub{parIndex,1}.metaData.meanHeadPitchList;
    plot(   plotThis(plotThis~=0)   ,'LineWidth',2)
    clear plotThis
end
ylim([-5 5])
ylabel('Head Pitch - PreTrial - Degrees on Sphere')
xlabel('pretrial gaze scene #')
legend(char(dataStruct.subjectName'))
tidyFigs
title('Mean AVG Head Pitch, last 250 samples of pre-trial scene');

%%plot  mean raw gaze pitch
figure('DefaultAxesFontSize',13)
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = dataStruct.sub{parIndex,1}.metaData.meanGazePitchList;
    plot(   plotThis(plotThis~=0)   ,'LineWidth',2)
    clear plotThis
end
ylim([-10 10])
ylabel('raw gaze Pitch - PreTrial - Degrees on Sphere')
xlabel('pretrial gaze scene #')
legend(char(dataStruct.subjectName'))
tidyFigs
title('Mean raw Gaze Pitch, last 250 samples of pre-trial scene');

%%plot  mean raw gaze Yaw
figure('DefaultAxesFontSize',13)
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = dataStruct.sub{parIndex,1}.metaData.meanGazeYawList;
    plot(   plotThis(plotThis~=0)   ,'LineWidth',2)
    clear plotThis
end
ylim([-10 10])
ylabel('raw gaze Yaw - PreTrial - Degrees on Sphere')
xlabel('pretrial gaze scene #')
legend(char(dataStruct.subjectName'))
tidyFigs
title('Mean raw Gaze Yaw, last 250 samples of pre-trial scene');




%%plot the number of scenes excluded by pretrial
figure('DefaultAxesFontSize',13)
plotThis = [];
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = [plotThis dataStruct.sub{parIndex,1}.metaData.excludedPretrialCount];
end
bar(plotThis)
ylim([0 40])
ylabel('# of scenes excluded by pretrial-fixation')
xlabel('participant')
xticks(1:length(params.subjectNames))
xtickangle(90)
xticklabels({params.subjectNames})
%legend(char(dataStruct.subjectName'))
tidyFigs
title('Count: scenes excluded by pretrial-fixation');



%%plot the number of scenes excluded by not scanned min thresh
figure('DefaultAxesFontSize',13)
plotThis = [];
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = [plotThis dataStruct.sub{parIndex,1}.metaData.excludedScannedCount];
end
bar(plotThis)
ylim([0 40])
ylabel('# of scenes excluded by eye-tracker-fail')
xlabel('participant')
xticks(1:length(params.subjectNames))
xtickangle(90)
xticklabels({params.subjectNames})
%legend(char(dataStruct.subjectName'))
tidyFigs
title('Count: scenes excluded by not scanning 2/3 scene');


%%plot the number of scenes excluded by low confidence
figure('DefaultAxesFontSize',13)
plotThis = [];
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    plotThis = [plotThis dataStruct.sub{parIndex,1}.metaData.lowConfidenceCount];
end
bar(plotThis)
ylim([0 40])
ylabel('# of scenes excluded by eye-tracker-fail')
xlabel('participant')
xticks(1:length(params.subjectNames))
xtickangle(90)
xticklabels({params.subjectNames})
%legend(char(dataStruct.subjectName'))
tidyFigs
title('Count: scenes excluded by low confidence');







%%plot the number of scenes excluded by low confidence | notScan | Drift
figure('DefaultAxesFontSize',13)
plotThis = [];
for parIndex = 1:size(params.subjectNames,1) % cyle through participants
    hold on
    value = dataStruct.sub{parIndex,1}.metaData.lowConfidenceCount + dataStruct.sub{parIndex,1}.metaData.excludedScannedCount + dataStruct.sub{parIndex,1}.metaData.excludedPretrialCount;
    plotThis = [plotThis value];
    clear value
end
bar(plotThis)
ylim([0 40])
ylabel('# of scenes excluded by eye-tracker-fail')
xlabel('participant')
xticks(1:length(params.subjectNames))
xtickangle(90)
xticklabels({params.subjectNames}) % TODO: check if this works in Octave
%legend(char(dataStruct.subjectName'))
tidyFigs
title('Count: scenes excluded by Confidence, Not Scanning, Pretrial-drift');



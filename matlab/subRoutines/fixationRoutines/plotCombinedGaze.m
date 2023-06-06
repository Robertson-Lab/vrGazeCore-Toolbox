function [] = plotCombinedGaze(resizeIM, yaw, pitch, subjectName, currentSceneText, confidence, paths, yaw_new, pitch_new, mean_yaw, mean_pitch, begin_fix, length_fix, fileID, params)
%plot combined head, raw, fix data on one image
%head in white, gaze black, fix cross

% make sure points are wrapped
[yaw,pitch] = wrapPointsEquirect( yaw, pitch, 360, 180, 0); 
[yaw_new,pitch_new] = wrapPointsEquirect( yaw_new, pitch_new, 360, 180, 0); 
[mean_yaw,mean_pitch] = wrapPointsEquirect( mean_yaw, mean_pitch, 360, 180, 0); 

% convert from degrees to nearest pixel
[headx, heady] = degreesToPixels(yaw,pitch,2000,1000);
[eyex, eyey] = degreesToPixels(yaw_new,pitch_new,2000,1000);
[fixx, fixy] = degreesToPixels(mean_yaw,mean_pitch,2000,1000);

%% Setup Fig
fig = figure('visible','on'), imshow(resizeIM);%start a new figure
title([subjectName currentSceneText 'rawHMDCenter']);

hold on
axis on;

%% Plot lines between head, eye
for i = 1:length(eyex)
    if abs(headx(i)-eyex(i)) < 250
        line([headx(i),eyex(i)],[heady(i),eyey(i)] , 'Color','white','LineWidth',0.01),
    end
end

%% Plot Head Dir by TIME

scatter(headx,heady,20,'w','filled');
scatter(headx,heady,10,'k','filled');
hold on
c = linspace(1,10,length(headx));
colormap('jet');
scatter(headx,heady,5,c,'filled');
%% Plot Eye Gaze by TIME
axis on;
scatter(eyex,eyey,50,'k','filled');
c = linspace(1,10,length(eyex));
scatter(eyex,eyey,20,c,'filled');

%% save temp plot
filename = [paths.subjectFixPlotsDir subjectName '_' currentSceneText '_combinedPlotTemp'];
outputFilename = [filename, '.jpg'];
f=getframe;
imwrite(f.cdata,outputFilename);
close
%% Plot Fixations by Duration, Spread

tempIM = imread(outputFilename);
fig = figure('visible',params.plotVisibility), imshow(tempIM); %pic_width x pic_width x 3
title([subjectName currentSceneText 'rawHMDCenter']);
hold on
delete(outputFilename)
axis on;

% Calculate Spread/Ecc
point_spread = pointSpread(yaw_new,pitch_new,mean_yaw,mean_pitch,begin_fix,length_fix,fileID);
spread_scaled = 80+3000*(point_spread - min(point_spread))  /  (max(point_spread) - min(point_spread));

% Plot fixations
scatter(fixx,fixy,75,'r','filled');
scatter(fixx,fixy,6.5*spread_scaled,'r+','LineWidth',2.5);
scatter(fixx,fixy,6*spread_scaled,length_fix,'+','LineWidth',1);
colormap('hot');
scatter(fixx,fixy,65,length_fix,'filled');
ax=gca;
ax.Visible = 'off';

%save combined
try
filename = [paths.subjectFixPlotsDir subjectName '_' currentSceneText '_combinedPlot'];
outputFilename = [filename, '.jpg'];
f=getframe;
imwrite(f.cdata,outputFilename);
catch 
end

end
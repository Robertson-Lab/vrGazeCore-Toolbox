function [point_spread] = plotFix(resizeIM, mean_yaw, mean_pitch, yaw_new, pitch_new, begin_fix, length_fix, subjectName, currentSceneText, fileID, duration, paths, params)
%plot raw fixation points coded by color onto each image

[fixx, fixy] = degreesToPixels(mean_yaw,mean_pitch,2000,1000);
        
% Calculate Spread/Ecc

point_spread = pointSpread(yaw_new,pitch_new,mean_yaw,mean_pitch,begin_fix,length_fix,fileID);

spread_scaled = 80+3000*(point_spread - min(point_spread))  /  (max(point_spread) - min(point_spread));

%% Plot by Duration
title1 = [subjectName currentSceneText 'fixation Duration'];

fig = figure('visible',params.plotVisibility), imshow(resizeIM); %pic_width x pic_width x 3
title(title1);
hold on;
axis on;

scatter(fixx,fixy,75,'r','filled');
scatter(fixx,fixy,6.5*spread_scaled,'r+','LineWidth',2.5);
scatter(fixx,fixy,6*spread_scaled,length_fix,'+','LineWidth',1);
colormap('hot');
scatter(fixx,fixy,80,duration,'filled');
ax=gca;
ax.Visible = 'off';
hold on;

filename = [paths.subjectFixPlotsDir subjectName '_' currentSceneText 'b_fixations_Duration'];
outputFilename = [filename, '.jpg'];
f=getframe;
imwrite(f.cdata,outputFilename);
hold off
end
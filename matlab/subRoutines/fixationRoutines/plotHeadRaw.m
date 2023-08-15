function [] = plotHeadRaw(resizeIM, yaw, pitch, subjectName, currentSceneText, confidence, paths, params)
%plot raw fixation points coded by color onto each image

plotConfidence = 0;

for i = 1:size(pitch)
        if pitch(i) > 90
            pitch(i,1) = pitch(i)-180;
        elseif pitch(i) < -90
            pitch(i,1) = pitch(i)+180;
        end
        if yaw(i) > 180
            yaw(i,1) = yaw(i)-360;
        elseif yaw(i) < -180
            yaw(i,1) = yaw(i)+360;
        end
end
[x, y] = degreesToPixels(yaw,pitch,2000,1000);
%% Plot by TIME

fig = figure('visible',params.plotVisibility), imshow(resizeIM); %pic_width x pic_width x 3
title([subjectName currentSceneText 'rawHMDCenter']);
hold on
axis on;
scatter(x,y,30,'w','filled');
scatter(x,y,12,'k','filled');
hold on
c = linspace(1,10,length(x));
colormap('jet');
scatter(x,y,7,c,'filled');

filename = [paths.subjectFixPlotsDir subjectName '_' currentSceneText 'a_rawHMDCenter'];
outputFilename = [filename, '.jpg'];
f=getframe;
imwrite(f.cdata,outputFilename);
hold off
close

%% Plot by Confidence
if plotConfidence == 1
    confidenceScaled = confidence;
    confidenceScaled(end) = 0;
    confidenceScaled(end-1) = 1;

    fig = figure('visible',params.plotVisibility), imshow(resizeIM); %pic_width x pic_width x 3
    title([subjectName currentSceneText 'rawConfidence']);
    hold on
    axis on;
    scatter(x,y,35,'k','filled');
    hold on
    colormap('winter');
    scatter(x,y,7,confidenceScaled,'filled');

    filename = [paths.subjectFixPlotsDir subjectName '_' currentSceneText 'a_rawConfidence'];
    outputFilename = [filename, '.jpg'];
    f=getframe;
    imwrite(f.cdata,outputFilename);
    hold off
    close
end
clear x y
end

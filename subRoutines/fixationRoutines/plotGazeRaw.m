function [] = plotGazeRaw(resizeIM, yaw_new, pitch_new, subjectName, currentSceneText, confidence, paths, params)
%plot raw fixation points coded by color onto each image
%JSM - 6/6/18

plotConfidence = 0;
[x, y] = degreesToPixels(yaw_new,pitch_new,2000,1000);

%% Plot by TIME
%figure%start a new figure
%title([subjectName currentSceneText 'raw']);

fig = figure('visible',params.plotVisibility), imshow(resizeIM); %pic_width x pic_width x 3
title([subjectName currentSceneText 'raw']);
hold on
axis on;
scatter(x,y,35,'k','filled');
hold on
c = linspace(1,10,length(x));
colormap('jet');
scatter(x,y,7,c,'filled');

filename = [paths.subjectFixPlotsDir subjectName '_' currentSceneText 'a_raw'];
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
    
    %figure%start a new figure for plotting by confidence
    %title([paths.subjectFixPlotsDir currentSceneText 'rawConfidence']);

    fig = figure('visible',params.plotVisibility), imshow(resizeIM); %pic_width x pic_width x 3
    title([paths.subjectFixPlotsDir currentSceneText 'rawConfidence']);
    hold on
    axis on;
    scatter(x,y,35,'k','filled');
    hold on
    %c = linspace(1,10,length(x));
    colormap('winter');
    scatter(x,y,7,confidenceScaled,'filled');

    filename = [subjectFixPlotsDir subjectName '_' currentSceneText 'a_rawConfidence'];
    outputFilename = [filename, '.jpg'];
    f=getframe;
    imwrite(f.cdata,outputFilename);
    hold off
    close
end

clear x y
end

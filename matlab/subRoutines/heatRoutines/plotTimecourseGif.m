%%script to plot gifs to easily compare
function plotTimecourseGif(inputSubjects, paths, params)
    subjectName = inputSubjects;

for sceneIdx = 1:length(params.sceneList)
    sceneName = params.sceneList(sceneIdx);
    timecoursePath = [ paths.projectTimecoursePlotsDir subjectName '/' char(sceneName) '/'];
    savePath = [paths.projectTimecoursePlotsDir subjectName '/'];

    timeImages = dir([timecoursePath '*.jpeg']);
    nImages = length(timeImages);

    savename = [savePath subjectName '_' char(sceneName) '_timecourse.gif'];

    h =  figure;

    for idx = 1:nImages
        if isfile([timecoursePath timeImages(idx).name]) == 1
            im2show = imread([timecoursePath timeImages(idx).name]);
            imshow(im2show)
    
            text(420,50,{[ subjectName ' scene:' char(sceneName) ' segment 01' ]},'vert','bottom','horiz','center','FontSize',15,'Color','blue');
            
            %make borderless
            set(gca,'units','pixels'); % set the axes units to pixels
            x = get(gca,'position'); % get the position of the axes
            set(gcf,'units','pixels'); % set the figure units to pixels
            y = get(gcf,'position'); % get the figure position
            %set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
            set(gcf,'position',[232 343 1108 565])
            set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
    
    
            frame = getframe(h); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
    
            if idx == 1
                imwrite(imind,cm,savename,'gif', 'DelayTime', params.delayTime, 'Loopcount',inf);
            else
            imwrite(imind,cm,savename,'gif', 'DelayTime', params.delayTime, 'WriteMode','append'); 
            end
        else
            fprintf('Scene heatmap %s not found, skipping to next scene heatmap.\n', [char(sceneName) '_' num2str(idx)])
        end
    end
    close all
    
    %remove .jpeg images
    if params.deleteTimecourseJPEG == 1
        rmdir (timecoursePath, 's')     
    end
end

end %end fn
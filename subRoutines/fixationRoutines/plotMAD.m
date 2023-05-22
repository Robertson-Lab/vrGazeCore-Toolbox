function [] = plotMAD(gazeX_deg,gazeY_deg,vel,madVel,time,swTime,endSWIdx,subjectName,currentSceneText,paths,params)
    % function to plot gaze position, velocity, and MAD velocity over time of trial  
    
    normTime = time - min(time);
    normSWTime = swTime - min(time);

    f = figure('Position',get(0,'Screensize'),'Visible',params.plotVisibility);
    
    tiledlayout(4,1)

    nexttile
    plot(normTime,gazeX_deg,"LineWidth",1)
    xlabel('Trial Time')
    ylabel('X Gaze Position (Viewport Deg)')
    set(gca, 'box', 'off')
    set(gca,'TickDir','out')
    title1 = [subjectName currentSceneText];
    title(title1)

    nexttile
    plot(normTime,gazeY_deg,"LineWidth",1)
    xlabel('Trial Time')
    ylabel('Y Gaze Position (Viewport Deg)')
    set(gca, 'box', 'off')
    set(gca,'TickDir','out')

    nexttile
    plot(normSWTime,vel,"LineWidth",1)
    xlabel('Trial Time')
    ylabel('Velocity (degs/sec)')
    set(gca, 'box', 'off')
    set(gca,'TickDir','out')

    nexttile
    plot(normSWTime(1:endSWIdx),madVel,"LineWidth",1)
    hold on
    plot(normSWTime(1:endSWIdx),repmat(params.minMad,length(swTime(1:endSWIdx))),'--','LineWidth',2,'Color','r')
    xlabel('Trial Time')
    ylabel('Mean Absolute Deviation (deg/sec)')
    set(gca, 'box', 'off')
    set(gca,'TickDir','out')

    filename = [paths.subjectMADPlotsDir subjectName '_' currentSceneText 'MAD-plots'];
    outputFilename = [filename, '.jpg'];
    saveas(gca,outputFilename);

    close
end
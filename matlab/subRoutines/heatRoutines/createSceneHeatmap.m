function [sceneHeat] = createSceneHeatmap(fixationEquiX, fixationEquiY, fixationDurations, trimFactor)

    sceneHeat = zeros(2000, 1000); % start a heat image where everywhere is zero (no heat)

    for i = 1:length(fixationDurations) % loop through fixations and mark each on heat image based on its duration weighting
        if isnan(fixationDurations(i))
            continue % at one point, the sceneDur length for a person was NaN - not sure why - problem stems from calculateFixations?
        end

        if fixationEquiY(i) < 1000*trimFactor+1 || fixationEquiY(i) > (1000-(1000*trimFactor))
            continue % if y coordinate is above or below trim region, skip.
        end
        if round(fixationEquiX(i)) == 0
            fixationEquiX(i) = 1; % no zero indexing in matlab
        end
        if round(fixationEquiY(i)) == 0
            fixationEquiY(i) = 1; % no zero indexing in matlab
        end

        % TODO: DOCUMENT THIS STEP
        sceneHeat(round(fixationEquiX(i)),round(fixationEquiY(i))) = sceneHeat(round(fixationEquiX(i)),round(fixationEquiY(i)))+fixationDurations(i);
    end

    warning('off', 'signal:check_order:InvalidOrderRounding') %turns off warning that pops up when running gaussianFilterEquirect
    sceneHeat = flipud(rot90(sceneHeat));
    sceneHeat = gaussianFilterEquirect(sceneHeat,200);

end
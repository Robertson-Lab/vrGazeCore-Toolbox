function [fixationDurationsFilteredScaled] = scaleFixationDurations(fixationDurations, params)

    % find the upper and lower bounds (top 95th percentile,
    % bottom 0.1 percentile
    if params.boundFiltering == 1
        upperBound = prctile(fixationDurations,95); % set upper bound to 95 percentile
        lowerBound = prctile(fixationDurations,0.1); % set lower bound to 0.1 percentile
    elseif params.boundFiltering == 0  % don't mess with bounds
        upperBound = prctile(fixationDurations,100);
        lowerBound = prctile(fixationDurations,0);
    end

    % based off of Henderson & Hayes (2018)
    % filter the durations based on the bounds to create duration-weighted maps 
    fixationDurationsFiltered = fixationDurations;
    fixationDurationsFiltered(fixationDurationsFiltered>upperBound) = upperBound;% set everything above upperbound to upperbound
    fixationDurationsFiltered(fixationDurationsFiltered<lowerBound) = lowerBound;% set everything below lowerbound to lowerbound

    % this normalizes fixation durations
    fixationDurationsFilteredScaled = (fixationDurationsFiltered - lowerBound) / (upperBound - lowerBound); % scale to 0 to 1
    fixationDurationsFilteredScaled=0.1+(fixationDurationsFilteredScaled*0.9); % scale to 0.1 to 1
end
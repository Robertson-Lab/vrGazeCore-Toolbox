%Function to load and resize all heat maps

function [mapOut] = trimMap(mapIn, trimFactor)
    [mapH,mapW] = size(mapIn);
    
    mapOut = mapIn((mapH*trimFactor+1):(mapH-(mapH*trimFactor)),:);

end

%% DP's Notes:
%  PURPOSE: wraps points around an equirectangluar projection
%  INPUTS: points in x-dim, y-dim, dimensions of equirect projections in x-dim, dimentions of equirect projections in y-dim, roundpoints
%  OUTPUTS: same variable names, just with equirect corrected

function [ inputX, inputY ] = wrapPointsEquirect( inputX, inputY, maxX, maxY, roundPoints)
%WRAPPOINTSEQUIRECT wrap points around an equirectangular projection
% Some fixation locations are off the screen, shift them so that
% they are in the correct position on the equirect projection.
% ~~~~~~~~~~~~~~~~~~
% Inputs:
% inputX: list of points in x-dimension
% inputY: list of points in y-dimension
% maxX: dimensions of equirect projection in x dimension
    % if equirect, it is imgW pixels
    % if degrees, it is 360 (yaw)
% maxY: dimensions of equirect projection in y dimension
    % if equirect, it is imgH pixels
    % if degrees, it is 180 (yaw)
% roundPoints: 1 if rounding (do this for pixel only for plotting purposes!
% ~~~~~~~~~~~~~~~~~~
if roundPoints == 1 
    inputY = round(inputY); % round to nearest pixel location
    inputX = round(inputX); % same for X coordinates
end

% if fixation is below screen, shift it 180d (1000px) in x dimension to wrap around pole
inputX(find(inputY<0 & inputX>(maxX/2))) = inputX(find(inputY<0 & inputX>(maxX/2)))-(maxX/2); 
inputX(find(inputY<0 & inputX<(maxX/2))) = inputX(find(inputY<0 & inputX<(maxX/2)))+(maxX/2); % if fixation is below screen, shift it back up from the bottom
% if fixation is below screen, shift it back up from the bottom
inputY(find(inputY<0)) = -inputY(find(inputY<0)); %if less than zero, make it positive y dimension

% if fixation is above screen, shift it 180d (1000px) in x dimension to wrap around pole
inputX(find(inputY>maxY & inputX>(maxX/2))) = inputX(find(inputY>maxY & inputX>(maxX/2)))-(maxX/2); 
inputX(find(inputY>maxY & inputX<(maxX/2))) = inputX(find(inputY>maxY & inputX<(maxX/2)))+(maxX/2); % if fixation is below screen, shift it back up from the bottom
% if fixation is above screen, shift it back down from the top
inputY(find(inputY>maxY)) = (maxY*2) - inputY(find(inputY>maxY)); %if less than zero, make it positive y dimension

% if fixation is to the left of screen, wrap around to the right + viceversa
inputX(find(inputX<0)) = inputX(find(inputX<0))+maxX;
inputX(find(inputX>maxX)) = inputX(find(inputX>maxX))-maxX;

end
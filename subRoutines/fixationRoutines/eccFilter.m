%% DP's Notes:
%  PURPOSE: removes eccentric points beyond given degree 
%  INPUTS: filter values
%  OUTPUTS: indices with invalid eccentricities

function [ds_out,ecc,percentRemoved] = eccFilter(ds_in,filterX,filterY,fovX,fovY,fileID)
%removes all eccentric points beyond given degrees (we think these errors in eyetracker reporting)
%   inputs:
%   ds - dataset
%   filterX,filterY - threshold 

    filterX = filterX/2; % go from 100 to +/- 50, since 0,0 is in center of screen
    filterY = filterY/2;
    pre_gazeX = ds_in(:,5);
    pre_gazeY = ds_in(:,6);
    pre_gazeXdeg = (pre_gazeX-0.5)*fovX; %convert X to degrees
    pre_gazeYdeg = (pre_gazeY-0.5)*fovY; %convert Y to degrees
    goodX = find(pre_gazeXdeg<filterX & pre_gazeXdeg>(filterX * -1));% find indices within X threshold
    goodY = find(pre_gazeYdeg<filterY & pre_gazeYdeg>(filterY * -1));% find indices within Y threshold
    goodXY = intersect(goodX, goodY);%find intersection

%     %jeff old method
%     ds_out = ds_in(goodXY,:);%apply filter to ds
%     
    %tlb 8-29-2019 - hacky way to find missing idxs
    invalidIdxs = find(~ismember(goodXY, [1:length(ds_in)]'));
    ds_in(invalidIdxs,end) = 1;%update invalid idx column
    ds_out = ds_in;
    
%Log how much was removed in command window
dsL = length(ds_in(:,7));
ecc = length(ds_in(:,7))-length(goodXY);
percentRemoved = 100 * ecc/dsL;
fprintf(fileID,'ECCENTRICITY FILTER: removed %d out of %d values (%.2f %%)\n', ecc, dsL, ecc/dsL);
fprintf('ECCENTRICITY FILTER: removed %d out of %d values (%.2f %%)\n', ecc, dsL, ecc/dsL);



end

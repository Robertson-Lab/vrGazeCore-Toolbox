function [ds_out,confidenceFlag,hiC,percentRemoved] = confidenceFilter(data_in,params,fileID)
%   removes all values where confidence is below threshold
%   inputs: ds (seventh column is confidence with our data)
%   cutoff threshold [decimal number 0-1] a good value is 0.6 as of pupil
%   v. 1.2. Newer versions seem to have lower confidence values so this can
%   be lower

confidence = data_in(:,7); %running confidence values for each timepoint
hiConf = find(confidence>params.minConfThresh); %find rows where confidence was higher than our confidence threshold

%dding preservation of low confidence indices to be used
%in fixation calculation
idxMat = zeros(size(data_in,1),1);

%find invalid rows and label
idxMat(find(confidence<=params.minConfThresh)) = 1;
ds_out = [data_in, idxMat];

%Log how much was removed in command window
dsL = length(data_in(:,7));
hiC = length(hiConf);

percentRemoved = (100*(dsL-hiC)/dsL);

if percentRemoved > params.maxConfPercent
    confidenceFlag = 1;
else 
    confidenceFlag = 0;
end

fprintf(fileID,'CONFIDENCE FILTER: removed %d out of %d values (%.2f %%)\n', dsL-hiC, dsL, percentRemoved);
fprintf('CONFIDENCE FILTER: removed %d out of %d values (%.2f %%)\n', dsL-hiC, dsL, percentRemoved);

end

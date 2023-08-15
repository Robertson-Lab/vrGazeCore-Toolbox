function [point_spread] = pointSpread(yaw_new,pitch_new,mean_yaw,mean_pitch,begin_fix,length_fix,fileID)
%see spread of points that make up each fixation
%mean arc Distance in degrees from center point

%   inputs:
%   a list of fixations
% a list of data

%out: 'spread' of each fixation mean arc distance from centroid

yaw_new = yaw_new - 180;
pitch_new = pitch_new - 90;
mean_yaw = mean_yaw - 180;
mean_pitch = mean_pitch - 90;

for i = 1:length(begin_fix)
    
    x_raws = yaw_new(begin_fix(i):begin_fix(i)+length_fix(i)-1);
    y_raws = pitch_new(begin_fix(i):begin_fix(i)+length_fix(i)-1);
    x_fix = mean_yaw(i);
    y_fix = mean_pitch(i);
    
    for d = 1:length(x_raws)

        % TLB - 6/4/23 --> changing argument order to correct order
        [distances(d),az] = distance(y_fix, x_fix, y_raws(d), x_raws(d));
%         [distances(d),az] = distance(x_fix,y_fix,x_raws(d),y_raws(d));
        
    end
    
    %tlb testing converting to pixels first for plotting reasons
    
    dist_pitch = mean(distances);
    dist_yaw = abs(mean(distances) * (1/cosd(y_fix))); % scale yaw by pitch

    [spread_yaw(i), spread_pitch(i)] = degreesToPixels(dist_yaw, dist_pitch, 2000, 1000); %then convert the distance of spread into pix
    
    
    point_spread(i) = mean(distances);
    
    
    clear x_raws y_raws distances
end




end

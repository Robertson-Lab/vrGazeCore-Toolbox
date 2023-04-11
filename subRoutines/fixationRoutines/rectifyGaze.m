%% DP's Notes:
%  PURPOSE: rectifies gaze data with head data
%  INPUTS: filtered gaze data, pitch, yaw
%  OUTPUTS: pitch_new, yaw_new

function [yaw_new,pitch_new] = rectifyGaze(gazeX_deg, roll_rad, gazeY_deg, pitch, yaw)
% JSM 12/13/18 -adding python fn that does this accurately
% Function to rotate all gaze points based on roll, rectifying gaze data with head data

pitch = pitch-90;
yaw = yaw-180;

for i = 1:size(gazeX_deg)
    % Calculate change in pitch and yaw accounting for rotation
    dx(i,1) = (gazeX_deg(i))*cos(roll_rad(i))-(gazeY_deg(i))*sin(roll_rad(i));
    dy(i,1) = (gazeY_deg(i))*cos(roll_rad(i))+(gazeX_deg(i))*sin(roll_rad(i));
    % shift pitch
    pitch_new(i,1) = pitch(i)-dy(i); %dy is subtracted since up is negative
    % scale yaw by pitch (because operating in latlon)
    dx(i) = dx(i) * cosd(pitch(i,1));
    yaw_new(i,1) = yaw(i)+dx(i);
end

usePython = 0; % don't use this for now but it is the right way accounting for gnomonic projection......
if usePython == 1
    if count(py.sys.path,'') == 0
        insert(py.sys.path,int32(0),'');
    end
    py.importlib.import_module('numpy');
    py.gnomonic2latlon.gnomonic2lat_lon()
    py.numpy.array([1,1])
    
end


pitch_new = pitch_new+90;
yaw_new = yaw_new+180;

end

function [time,pitch,yaw,roll_rad,gazeX_deg,gazeY_deg,gaze_ecc,confidence] = parseDS(ds_filtered,params)
%Function to pull out and name useful variables from ds_filtered
%Also calculates "gaze eccentricity" (distance from 0,0)
%
% NOTE - In this context, all variables still refer to eye data (not
% eye+head)

	% Name variables
	time = ds_filtered(:,1); % experiment clock time from Unity
	pitch = ds_filtered(:,2); % 'x' in unity is pitch and negative when looking upwards
	yaw = ds_filtered(:,3);
	roll = ds_filtered(:,4);
    
    if params.headsetType == 3
        roll_rad = [];
        gazeX_deg = [];
        gazeY_deg = [];
        gaze_ecc = [];
        confidence = [];
        return;
        
    end
    
	gazeX = ds_filtered(:,5);
	gazeY = ds_filtered(:,6);
    confidence = ds_filtered(:,7);
	% Convert roll into radians
	roll_rad = deg2rad(roll);
    
    % shift to positive domain for easier scaling
    pitch = pitch+90; % shift from -90to90 to 0to180
    yaw = yaw+180; % shift from -180to180 to 0to360

    %  Ehringer (2019): we converted the x (and y) gaze points of the raw samples from screen coordinates in pixels, 
    %  to spherical angles in degree (with a reference system centered on the subject):
    %
    %           Bx = 2 � atan2(px � m, d)
    %
    %  where Bx denotes the azimuth angle (equivalent to the horizontal position) of the gaze points 
    %  in visual degrees from the monitor center
    %
    %       - px: horizontal position relative to the center of the monitor in pixel, 
    %       - m: unit conversion of pixel to mm of the monitor, 
    %       - d: the distance to the monitor in mm.
    
    % X,Y coordinates are presented on a plane ~0.8 units from the
    % camera
    
    if params.headsetType == 0 % DK2
        
        % normalize center of the screen (eye data) to 0,0. because y
        % calibration coord is slightly below center of screen, shift a
        % bit more than X.
        
        gazeX = gazeX-0.5;
        gazeY = gazeY-0.505;
        
        gazeX_deg = rad2deg ( (2*params.fovX/params.maxFOV) * atan2(  gazeX , 0.8) );
        gazeY_deg = rad2deg ( (2*params.fovY/params.maxFOV) * atan2(  gazeY , 0.8));

    elseif params.headsetType == 2
        if params.gazeType == 1 % Vive Eye - 2D tracking

            gazeY = gazeY+0.005;
            gazeX = params.maxFOV/params.fovX * asin(-1*gazeX); %though using tan works too
            gazeY = params.maxFOV/params.fovY * asin(gazeY);
            gazeX_deg = rad2deg(gazeX);
            gazeY_deg = rad2deg(gazeY);
            
        else % Vive Eye - 3D tracking
            
            gazeX = gazeX-0.5;
            gazeY = gazeY-0.5;
            gazeX_deg = rad2deg ( (2*params.fovX/params.maxFOV) * atan2(gazeX, 0.55));
            gazeY_deg = rad2deg ( (2*params.fovY/params.maxFOV) * atan2(gazeY, 0.55));
            
        end

    end
    
	% CALCULATE GAZE ECCENTRICITY (D) [Distance from 0,0 (degrees)]
    [gaze_ecc,~] = distance(0,0,gazeX_deg,gazeY_deg); %how eccentric gaze was from the viewport at a given point
end

%% DP's Notes:
%  PURPOSE: pulls out & names useful variables from parData (participant's raw data matrix); calculates gaxe eccentricity & drift correction
%  INPUTS: parData
%  OUTPUTS: 'raw...' variables 

function [rawSceneData] = parseRaw2struct(parData,fovX,fovY,preTrialScene,params,preshiftX,preshiftY)
%Function to pull out and name useful variables from parData
%Also calculates "gaze eccentricity" (distance from 0,0)
%also performs drift correction
	% Name variables
	rawSceneData = struct;
    
    rawTime = parData(:,1);
	rawHeadPitch = parData(:,2); % 'x' in unity is pitch and negative when looking upwards
	rawHeadYaw = parData(:,3);
	rawHeadRoll = parData(:,4);
    
    if params.headsetType == 3
        rawSceneData.rawTime = rawTime;
        rawSceneData.rawHeadPitch = rawHeadPitch;
        rawSceneData.rawHeadYaw = rawHeadYaw;
        rawSceneData.rawHeadRoll = rawHeadRoll;
        return;
    end
    
	rawViewportX_norm = parData(:,5);
	rawViewportY_norm = parData(:,6);
    rawConfidence = parData(:,7);
	% Convert roll into radians
	rawHeadRoll_rad = deg2rad(rawHeadRoll);
    
%     %shift to positive domain for easier scaling
%     pitch = pitch+90; % shift from -90to90 to 0to180
%     yaw = yaw+180; % shift from -180to180 to 0to360

    if params.unityProjectVersion == 1
        
        %  tlb 9-02-19: useful line explaining the conversion below
        %
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
        
%         % tlb - new method?

        % X,Y coordinates are presented on a plane ~0.8 units from the
        % camera
        
        if params.headsetType == 0 % DK2
            
            rawViewportX_norm = rawViewportX_norm-0.5;
            rawViewportX_deg = rad2deg ( (2*params.fovX/params.maxFOV) * atan2(  rawViewportX_norm , 0.8) );
            
            rawViewportY_norm = rawViewportY_norm-0.505;
            rawViewportY_deg = rad2deg ( (2*params.fovY/params.maxFOV) * atan2(  rawViewportY_norm , 0.8));

        elseif params.headsetType == 2
            

            
            if params.gazeType == 1
                rawViewportX_norm(rawViewportX_norm >= 0.25) = 1.15*rawViewportX_norm(rawViewportX_norm>=0.25);
                rawViewportX_norm(rawViewportX_norm <= -0.25) = 1.15*rawViewportX_norm(rawViewportX_norm<=-0.25);
                rawViewportY_norm(rawViewportY_norm >= 0.05) = 1.05*rawViewportY_norm(rawViewportY_norm>=0.05);
                rawViewportY_norm(rawViewportY_norm <= -0.05) = 1.05*rawViewportY_norm(rawViewportY_norm<=-0.05);

                rawViewportX_deg = rad2deg(rawViewportX_norm);
                rawViewportY_deg = rad2deg(rawViewportY_norm);
%                 gazeX_deg = rad2deg(-1.1*gazeX);
%                 gazeY_deg = rad2deg(gazeY);
                
            else
                rawViewportX_norm = rawViewportX_norm-0.5;
                rawViewportY_norm = rawViewportY_norm-0.5;
                rawViewportX_deg = rad2deg ( (2*params.fovX/params.maxFOV) * atan2(  rawViewportX_norm , 0.55) );

                rawViewportY_deg = rad2deg ( (2*params.fovY/params.maxFOV) * atan2(  rawViewportY_norm , 0.55));
                
            end

        end
    else
        
        % Use arctan(2x)
        rawViewportX_norm = rawViewportX_norm-0.5;
        rawViewportX_deg( rawViewportX_norm>=0 ,: ) = rad2deg ( 1.9 * atan( 1* rawViewportX_norm( rawViewportX_norm>=0 ) ) );
        rawViewportX_deg( rawViewportX_norm< 0 ,:) = - rad2deg ( 1.9 * atan( -1* rawViewportX_norm( rawViewportX_norm< 0 ) ) );

        rawViewportY_norm = rawViewportY_norm-0.505;
        rawViewportY_deg( rawViewportY_norm>=0 ,: ) = rad2deg ( 1.9 * atan( 1 * rawViewportY_norm( rawViewportY_norm>=0 ) ) );
        rawViewportY_deg( rawViewportY_norm< 0 ,:) = - rad2deg ( 1.9 * atan( -1* rawViewportY_norm( rawViewportY_norm< 0 ) ) );
    end
    
%     if params.unityProjectVersion == 1
%     % Use arctan(x) + constants to convert gaze from norm VPS to degrees
%         rawViewportX_norm = rawViewportX_norm-0.5;
%         rawViewportX_deg( rawViewportX_norm>=0 ,: ) = rad2deg ( 1.85 * atan( 1.3* rawViewportX_norm( rawViewportX_norm>=0 ) ) );
%         rawViewportX_deg( rawViewportX_norm< 0 ,:) = - rad2deg ( 1.85 * atan( -1.3* rawViewportX_norm( rawViewportX_norm< 0 ) ) );
% 
%         rawViewportY_norm = rawViewportY_norm-0.505;
%         rawViewportY_deg( rawViewportY_norm>=0 ,: ) = rad2deg ( 2.1 * atan( 1.3 * rawViewportY_norm( rawViewportY_norm>=0 ) ) );
%         rawViewportY_deg( rawViewportY_norm< 0 ,:) = - rad2deg ( 2.1 * atan( -1.3* rawViewportY_norm( rawViewportY_norm< 0 ) ) );
%     else
%         % Use arctan(2x)
%         rawViewportX_norm = rawViewportX_norm-0.5;
%         rawViewportX_deg( rawViewportX_norm>=0 ,: ) = rad2deg ( 1.9 * atan( 1* rawViewportX_norm( rawViewportX_norm>=0 ) ) );
%         rawViewportX_deg( rawViewportX_norm< 0 ,:) = - rad2deg ( 1.9 * atan( -1* rawViewportX_norm( rawViewportX_norm< 0 ) ) );
% 
%         rawViewportY_norm = rawViewportY_norm-0.505;
%         rawViewportY_deg( rawViewportY_norm>=0 ,: ) = rad2deg ( 1.9 * atan( 1 * rawViewportY_norm( rawViewportY_norm>=0 ) ) );
%         rawViewportY_deg( rawViewportY_norm< 0 ,:) = - rad2deg ( 1.9 * atan( -1* rawViewportY_norm( rawViewportY_norm< 0 ) ) );
%     end
%     

    if preTrialScene == 0 && params.driftCorrection == 1
        rawViewportX_deg_driftcorrect = -preshiftX + rawViewportX_deg ;
        rawViewportY_deg_driftcorrect = -preshiftY + rawViewportY_deg ;
    end

    rawSceneData.rawTime = rawTime;
    rawSceneData.rawHeadPitch = rawHeadPitch;
    rawSceneData.rawHeadYaw = rawHeadYaw;
    rawSceneData.rawHeadRoll = rawHeadRoll;
    rawSceneData.rawViewportX_norm = rawViewportX_norm;
    rawSceneData.rawViewportY_norm = rawViewportY_norm;
    rawSceneData.rawConfidence = rawConfidence;
    rawSceneData.rawHeadRoll_rad = rawHeadRoll_rad;
    rawSceneData.rawViewportX_deg = rawViewportX_deg;
    rawSceneData.rawViewportY_deg = rawViewportY_deg; 
    
    if preTrialScene == 0 && params.driftCorrection == 1
        rawSceneData.rawViewportX_deg_driftcorrect = rawViewportX_deg_driftcorrect;
        rawSceneData.rawViewportY_deg_driftcorrect = rawViewportY_deg_driftcorrect;
    end

end
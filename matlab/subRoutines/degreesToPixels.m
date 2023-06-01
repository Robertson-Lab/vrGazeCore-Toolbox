function [x_out,y_out] = degreesToPixels(yaw_in,pitch_in,imgW,imgH)
%scales degrees to coordinates in pixels on an equirectangular projection
%   inputs: yaw (degrees), pitch (degrees), image width (pixels), image
%   height (pixels)

% scaling from -90to90;-180to180 OLD
% x_out = ((yaw_in * imgW/2) / 180) + imgW/2;
% y_out = ((pitch_in * imgH/2) / 90) + imgH/2;

%scaling from 0to180;0to360 NEW
x_out = ((yaw_in * imgW) / 360);
y_out = ((pitch_in * imgH) / 180);

% round to nearest pixel
x_out = round(x_out);
y_out = round(y_out);

end
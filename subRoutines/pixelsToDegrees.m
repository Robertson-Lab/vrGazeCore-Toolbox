function [yaw_out,pitch_out] = pixelsToDegrees(x_in,y_in,imgW,imgH)
%scales degrees to coordinates in pixels on an equirectangular projection
%   outputs: yaw (degrees), pitch (degrees)
%   inputs: image width (pixels), image
%   height (pixels)

% scaling from -90to90;-180to180 OLD
% x_out = ((yaw_in * imgW/2) / 180) + imgW/2;
% y_out = ((pitch_in * imgH/2) / 90) + imgH/2;

%scaling from 0to180;0to360 NEW
yaw_out = (x_in *360) / imgW;
pitch_out = (y_in * 180) / imgH;

% % round to nearest pixel
% x_in = round(x_in);
% y_in = round(y_in);
end
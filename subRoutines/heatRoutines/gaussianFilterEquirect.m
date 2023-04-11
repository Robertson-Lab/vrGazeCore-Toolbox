%% gaussianFilterEquirect
% function to go from points(e.g., fixations) on an equirectangular image to gaussian filtered
% heatmap where the gaussian width is variable depending on the pitch to account for equirectangular distortion 
% Can also be applied to smoothing pixel maps?

% try doing rows before columns

function [imgOut] = gaussianFilterEquirect(imgIn,baseGaussWidth)
%% Test that this works with fake fixations
% baseGaussWidth = 200;
% % 
% imgIn = zeros(1000,2000);
% % for i = 1:20
% %     imgIn(randi(1000) ,randi(2000)) = 1;
% % end
% imgIn(500 ,1000) = .8;
% imgIn(800 ,250) = .4;
% imgIn(300 ,1300) = .3;
% imgIn(50 ,50) = .2;
% imgIn(999 ,1000) = .1;
% 
% imgIn(995 ,995) = .1;
% imgIn(990 ,980) = .1;
% imgIn(980 ,990) = .1;



% get image H,W
imageW = size(imgIn,2);
imageH = size(imgIn,1);



%% Pad & smooth the image - Horizontal steps: pad, variable width gaussian
% pad BOTH L/R here (no need to reflect for the horizontal step)
imgIn = [imgIn(:,1+(3*imageW/4):end) imgIn imgIn(:,1:(imageW/4))];


%make copy of img that is flipped
imgRowRev = fliplr(imgIn);


%loop through rows of original
for i = 1:imageH
    % get pitch in degrees from pixels
    [~,pitch_out] = pixelsToDegrees(i,i,imageW,imageH);
    pitchDegrees = pitch_out-90;
    % set variable size of gaussian width based on the pitch
    variableGaussWidth = baseGaussWidth*(1/cosd(pitchDegrees));
%     if variableGaussWidth == inf || variableGaussWidth>10000
%         variableGaussWidth = 10000;
%     end
    if variableGaussWidth == inf
        variableGaussWidth = 1000000;
    end
    w = gausswin(variableGaussWidth); %get gaussianwindow
    w = w(round(length(w)/2):end); % get half of the gaussian window
    imgIn(i,:) = filter(w,1,imgIn(i,:)); % filter each row
end

%loop again for the reversed
for i = 1:imageH
    % get pitch in degrees from pixels
    [~,pitch_out] = pixelsToDegrees(i,i,imageW,imageH);
    pitchDegrees = pitch_out-90;
    variableGaussWidth = baseGaussWidth*(1/cosd(pitchDegrees));

%     if variableGaussWidth == inf || variableGaussWidth>10000
%         variableGaussWidth = 10000;
%     end
    if variableGaussWidth == inf
        variableGaussWidth = 1000000;
    end
    w = gausswin(variableGaussWidth);
    w = w(round(length(w)/2):end);
    imgRowRev(i,:) = filter(w,1,imgRowRev(i,:));
end

imgRowRev = [ zeros(imageH,1) imgRowRev]; % since they overlap by one column, trim one of the columns
imgRowRev(:,end) = []; % add an extra to the other side so the matrices are the same size
imgRowRev = fliplr( imgRowRev); % flip it back


imgOut = imgIn + imgRowRev; % add them together

% unpad L/R here
imgOut = imgOut(:,1+(imageW/4):end-(imageW/4));






%% Pad & smooth the image - Vertical steps: pad, fixed width gaussian
% We create padding out of actual image data, but first we reflect it 
% imgInUD = flipud(imgOut);
% imgInUDLR = fliplr(imgInUD);


% Then we pad the top and bottom. 
% imgCol = [ imgInUDLR' imgOut' imgInUDLR'];
% imgCol = imgCol';
% 
% imgColUD = flipud(imgCol);
% 

% make gaussian window of fixed width
w = gausswin(baseGaussWidth);
w = w(round(length(w)/2):end); % get half of the gaussian window

imgOutUD = flipud(imgOut);
% loop through columns w/ fixed-width gaussian (twice, b/c asymmetrical)
for i = 1:imageW
    imgOut(:,i) = filter(w,1,imgOut(:,i));
end

% loop through the reverse image
for i = 1:imageW
    imgOutUD(:,i) = filter(w,1,imgOutUD(:,i));
end

imgOutUD = [ zeros(imageW,1) imgOutUD']; % since they overlap by one column, trim one of the columns
imgOutUD = imgOutUD';
imgOutUD(end,:) = []; % add an extra to the other side so the matrices are the same size
imgOutUD = flipud( imgOutUD); % flip it back


% unpad top and bottom here 
imgOut = imgOutUD + imgOut;
%imagesc(imgOut) %commented out to supress heatmap showing

% maximgout = max(max(imgOut(200:800,:)));
% 
% imgOut(imgOut>=maximgout) = maximgout;

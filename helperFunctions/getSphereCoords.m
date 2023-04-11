function [sampMat,indicesOut] = getSphereCoords(sampRate,imSize)

sphereSampling = SpiralSampleSphere(sampRate);
sphereSampling = rad2deg(sphereSampling);


sphereSampling(:,:) = sphereSampling(:,:) + repmat([180, 90], size(sphereSampling,1),1); %
sphereSampling(:,1) = (sphereSampling(:,1)/180)*imSize(2)/2;
sphereSampling(:,2) = (sphereSampling(:,2)/90)*imSize(1)/2;

indicesOut = round(sphereSampling);
indicesOut(indicesOut==0) = 1;

sampMat = zeros(imSize);
for i = 1:length(indicesOut)
    sampMat(indicesOut(i,2),indicesOut(i,1)) = 1;
end
sampMat = sampMat(:);
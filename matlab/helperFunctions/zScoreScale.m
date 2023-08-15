%Function to zscore each features and rescale them to -1 to 1
%9/15/18
%input: X - MxN array of all features (M = points, N = features)
%output: Xn - zscored and rescaled to -1 to 1 range OR 0 to 1
function [Xn] = zScoreScale(X)
[m,n] = size(X);% get size of input array
for i = 1:n
    X(~isnan(X(:,i)),i) = zscore(X(~isnan(X(:,i)),i)); %zscore
    X(:,i) = X(:,i) - min(X(:,i)); %rescale
    range = max(X(:,i)) - min(X(:,i));
    %X(:,i) = ( (X(:,i)/range)  *2 ) - 1; %rescale to range -1 to 1
    X(:,i) = (X(:,i)/range); %rescale to range 0 to 1
end
Xn = X;
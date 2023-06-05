function [ outputDVA ] = viewportNorm2dva(inputNorm )
%VIEWPORTNORM2DVA Summary of this function goes here
%   Detailed explanation goes here

coeffs = loadn2dcoef;

outputDVA = (coeffs(1)*power(inputNorm,3) + coeffs(2)*power(inputNorm,2) + coeffs(3)*power(inputNorm,1) + 0); %returns radius
end


function [ coeffs ] = loadn2dcoef(  )
%LOADN2DCOEF Summary of this function goes here
%   Detailed explanation goes here

norm2dva = csvread('norm2DVA.csv');
coeffs = polyfit(norm2dva(:,1), norm2dva(:,2), 3);

end


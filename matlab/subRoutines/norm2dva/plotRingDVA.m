% shitty script to plot table of ring values in Unity to DVA values
% and fit a function to it

clear all; close all;

%load the csv
norm2dva = csvread('norm2DVA.csv');

f = figure; %grab a fig
hold on
% myFit = fit(norm2dva(:,1), norm2dva(:,2), 'poly3'); % x = dva, y = ring size
% f, plot(myFit, norm2dva(:,1), norm2dva(:,2));
ylim([0 100])
xlim([0 0.5])





coeffs = polyfit(norm2dva(:,1), norm2dva(:,2), 3);


x=-2:.1:2;
p=[coeffs(1) coeffs(2) coeffs(3) 0] % polynomial function
plot(x,polyval(p,x))

y = rad2deg( atan(x*2) );
plot(x,y);

norm2dva_vive = norm2dva(1:end-1,:);
norm2dva_vive(:,1) = norm2dva_vive(:,1).*(80/90);

viveCoeffs = [];
viveCoeffs = polyfit(norm2dva_vive(:,1), norm2dva_vive(:,2), 3);
x_vive=-2:.1:2;
p_vive=[viveCoeffs(1) viveCoeffs(2) viveCoeffs(3) 0] % polynomial function
plot(x_vive,polyval(p_vive,x_vive))

%test the function
testDVA = 0.04; % diameter
(coeffs(1)*power(testDVA,3) + coeffs(2)*power(testDVA,2) + coeffs(3)*power(testDVA,1) + 0) %returns radius


